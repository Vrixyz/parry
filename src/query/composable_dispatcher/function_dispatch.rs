use core::any::TypeId;
use std::collections::HashMap;

use log::warn;

use crate::{
    math::{Isometry, Real},
    query::QueryDispatcher,
    shape::{RoundShape, Shape},
};

use super::ComposableQueryDispatcher;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DispatcherTypeKey(pub TypeId, pub TypeId);

pub trait IntersectionWithDispatcher: Send + Sync {
    fn intersection_with_dispatcher(
        &self,
        dispatcher: &ComposableQueryDispatcher,
        pos12: &Isometry<Real>,
        s1: &dyn Shape,
        s2: &dyn Shape,
    ) -> Result<bool, ()>;
}
pub struct BoxedIntersectionWithDispatcher(Box<dyn IntersectionWithDispatcher>);
impl IntersectionWithDispatcher for BoxedIntersectionWithDispatcher {
    fn intersection_with_dispatcher(
        &self,
        dispatcher: &ComposableQueryDispatcher,
        pos12: &Isometry<Real>,
        s1: &dyn Shape,
        s2: &dyn Shape,
    ) -> Result<bool, ()> {
        self.0
            .intersection_with_dispatcher(dispatcher, pos12, s1, s2)
    }
}
impl IntersectionWithDispatcher for fn(&Isometry<Real>, &dyn Shape, &dyn Shape) -> bool {
    fn intersection_with_dispatcher(
        &self,
        _dispatcher: &ComposableQueryDispatcher,
        pos12: &Isometry<Real>,
        s1: &dyn Shape,
        s2: &dyn Shape,
    ) -> Result<bool, ()> {
        Ok(self(pos12, s1, s2))
    }
}
/*
impl IntersectionWithDispatcher
    for fn(&ComposableQueryDispatcher, &Isometry<Real>, &dyn Shape, &dyn Shape) -> bool
{
    fn intersection_with_dispatcher(
        &self,
        dispatcher: &ComposableQueryDispatcher,
        pos12: &Isometry<Real>,
        other: &dyn Shape,
    ) -> bool {
        self(dispatcher, pos12, other, other)
    }
}*/

impl<F> IntersectionWithDispatcher for F
where
    F: Fn(&ComposableQueryDispatcher, &Isometry<Real>, &dyn Shape, &dyn Shape) -> Result<bool, ()>
        + Send
        + Sync,
{
    fn intersection_with_dispatcher(
        &self,
        dispatcher: &ComposableQueryDispatcher,
        pos12: &Isometry<Real>,
        s1: &dyn Shape,
        s2: &dyn Shape,
    ) -> Result<bool, ()> {
        self(dispatcher, pos12, s1, s2)
    }
}

pub struct FunctionDispatch {
    /// Map for the function to call when both parameters are known types.
    pub functions: HashMap<DispatcherTypeKey, Box<dyn IntersectionWithDispatcher + 'static>>,
}

impl core::fmt::Debug for FunctionDispatch {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FunctionDispatch")
            .field("known_both", &self.functions.keys())
            .finish()
    }
}

impl FunctionDispatch {
    pub fn new() -> Self {
        Self {
            functions: HashMap::new(),
        }
    }
    pub fn add_function_known_12<S1: Shape, S2: Shape>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &S1, &S2) -> bool + 'static + Send + Sync + Copy,
    ) {
        if self
            .functions
            .insert(
                DispatcherTypeKey(TypeId::of::<S1>(), TypeId::of::<S2>()),
                Box::new(
                    move |_dispatcher: &ComposableQueryDispatcher,
                          pose: &Isometry<Real>,
                          s1: &dyn Shape,
                          s2: &dyn Shape| {
                        let shape1 = s1.as_shape::<S1>().ok_or(())?;
                        let shape2 = s2.as_shape::<S2>().ok_or(())?;
                        Ok(inner(pose, shape1, shape2))
                    },
                ),
            )
            .is_some()
        {
            warn!(
                "Overwriting function for types {:?} and {:?}",
                TypeId::of::<S1>(),
                TypeId::of::<S2>()
            );
        }
    }

    pub fn add_function_known_1x<S1: Shape>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &S1, &dyn Shape) -> bool + 'static + Send + Sync + Copy,
        tids: Vec<TypeId>,
    ) {
        for type_s2 in tids {
            self.add_both_combinations(inner, TypeId::of::<S1>(), type_s2);
            self.add_both_combinations(inner, TypeId::of::<RoundShape<S1>>(), type_s2);
        }
    }

    pub fn add_function_dyn_dispatcher(
        &mut self,
        inner: impl Fn(&ComposableQueryDispatcher, &Isometry<Real>, &dyn Shape, &dyn Shape) -> bool
            + 'static
            + Send
            + Sync
            + Copy,
        type_s1: TypeId,
        type_s2: TypeId,
    ) {
        self.add_raw_function(
            type_s1,
            type_s2,
            move |dispatcher: &ComposableQueryDispatcher,
                  pose: &Isometry<Real>,
                  s1: &dyn Shape,
                  s2: &dyn Shape| { Ok(inner(dispatcher, pose, s1, s2)) },
        );

        self.add_raw_function(
            type_s2,
            type_s1,
            move |dispatcher: &ComposableQueryDispatcher,
                  pose: &Isometry<Real>,
                  s1: &dyn Shape,
                  s2: &dyn Shape| { Ok(inner(dispatcher, &pose.inverse(), s2, s1)) },
        );
    }

    fn add_both_combinations<S1: Shape>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &S1, &dyn Shape) -> bool + 'static + Send + Sync + Copy,
        type_s1: TypeId,
        type_s2: TypeId,
    ) {
        self.add_raw_function(
            type_s1,
            type_s2,
            move |_dispatcher: &ComposableQueryDispatcher,
                  pose: &Isometry<Real>,
                  s1: &dyn Shape,
                  s2: &dyn Shape| {
                let shape1 = s1.as_shape::<S1>().ok_or(())?;
                Ok(inner(pose, shape1, s2))
            },
        );
        // add the inverse function
        self.add_raw_function(
            type_s2,
            type_s1,
            move |_dispatcher: &ComposableQueryDispatcher,
                  pose: &Isometry<Real>,
                  s1: &dyn Shape,
                  s2: &dyn Shape| {
                let shape2 = s2.as_shape::<S1>().ok_or(())?;
                Ok(inner(&pose.inverse(), shape2, s1))
            },
        );
    }

    pub fn add_raw_function(
        &mut self,
        type_s1: TypeId,
        type_s2: TypeId,
        function: impl IntersectionWithDispatcher + 'static + Send + Sync + Copy,
    ) {
        if self
            .functions
            .insert(DispatcherTypeKey(type_s1, type_s2), Box::new(function))
            .is_some()
        {
            warn!(
                "Overwriting function for types {:?} and {:?}",
                type_s1, type_s2
            );
        }
    }

    fn add_both_combinations_dispatcher<S1: Shape>(
        &mut self,
        inner: impl Fn(&ComposableQueryDispatcher, &Isometry<Real>, &S1, &dyn Shape) -> bool
            + 'static
            + Send
            + Sync
            + Copy,
        type_s1: TypeId,
        type_s2: TypeId,
    ) {
        self.add_raw_function(
            type_s1,
            type_s2,
            move |dispatcher: &ComposableQueryDispatcher,
                  pose: &Isometry<Real>,
                  s1: &dyn Shape,
                  s2: &dyn Shape| {
                let shape1 = s1.as_shape::<S1>().ok_or(())?;
                Ok(inner(dispatcher, pose, shape1, s2))
            },
        );
        // add the inverse function
        self.add_raw_function(
            type_s2,
            type_s1,
            move |dispatcher: &ComposableQueryDispatcher,
                  pose: &Isometry<Real>,
                  s1: &dyn Shape,
                  s2: &dyn Shape| {
                let shape2 = s2.as_shape::<S1>().ok_or(())?;
                Ok(inner(dispatcher, &pose.inverse(), shape2, s1))
            },
        );
    }

    pub fn dispatch(
        &self,
        dispatcher: &ComposableQueryDispatcher,
        pose: &Isometry<Real>,
        s1: &dyn Shape,
        s2: &dyn Shape,
    ) -> Result<bool, ()> {
        let key = DispatcherTypeKey(s1.type_id(), s2.type_id());
        if let Some(func) = self.functions.get(&key) {
            return func.intersection_with_dispatcher(dispatcher, pose, s1, s2);
        }

        Err(())
    }
}

use core::any::TypeId;
use std::collections::HashMap;

use log::warn;

use crate::{
    math::{Isometry, Real},
    query::PointQuery,
    shape::{Shape, SimdCompositeShape, SupportMap, TypedSimdCompositeShape},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DispatcherTypeKey(pub TypeId, pub TypeId);

pub struct FunctionDispatch<'c, 'a: 'c, 'b: 'c> {
    /// Map for the function to call when both parameters are known types.
    pub known_both: HashMap<
        DispatcherTypeKey,
        Box<
            dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()>
                + Send
                + Sync
                + 'c,
        >,
    >,
    /// Map for the functions to try when first parameters is unknown.
    pub known_first: HashMap<
        TypeId,
        Vec<
            Box<
                dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()>
                    + Send
                    + Sync
                    + 'c,
            >,
        >,
    >,
    /// Map for the functions to try when second parameters is unknown.
    pub known_second: HashMap<
        TypeId,
        Vec<
            Box<
                dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()>
                    + Send
                    + Sync
                    + 'c,
            >,
        >,
    >,
    /// Map for the functions to try when both parameters are unknown.
    pub known_none: Vec<
        Box<
            dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()>
                + Send
                + Sync
                + 'c,
        >,
    >,
}

impl<'c, 'a: 'c, 'b: 'c> core::fmt::Debug for FunctionDispatch<'c, 'a, 'b> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FunctionDispatch")
            .field("known_both", &self.known_both.keys())
            .field("known_first", &self.known_first.keys())
            .field("known_second", &self.known_second.keys())
            .field("known_none", &self.known_none.len())
            .finish()
    }
}

impl<'c, 'a: 'c, 'b: 'c> FunctionDispatch<'c, 'a, 'b> {
    pub fn new() -> Self {
        Self {
            known_both: HashMap::new(),
            known_first: HashMap::new(),
            known_second: HashMap::new(),
            known_none: Vec::new(),
        }
    }

    pub fn add_function_known_12<S1: Shape + 'a, S2: Shape + 'b>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a S1, &'b S2) -> bool + 'c + Send + Sync + Copy,
    ) {
        if self
            .known_both
            .insert(
                DispatcherTypeKey(TypeId::of::<S1>(), TypeId::of::<S2>()),
                Box::new(
                    move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
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

    pub fn add_function_known_1<S1: Shape + 'a>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a S1, &'b dyn Shape) -> bool + 'c + Send + Sync + Copy,
    ) {
        let dispatch_second = self.known_first.entry(TypeId::of::<S1>()).or_insert(vec![]);
        dispatch_second.push(Box::new(
            move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                let shape1 = s1.as_shape::<S1>().ok_or(())?;
                let shape2 = s2;
                Ok(inner(pose, shape1, shape2))
            },
        ));
    }

    pub fn add_function_known_2<S2: Shape + 'b>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a dyn Shape, &'b S2) -> bool + 'c + Send + Sync + Copy,
    ) {
        let dispatch_second = self.known_first.entry(TypeId::of::<S2>()).or_insert(vec![]);
        dispatch_second.push(Box::new(
            move |pose: &Isometry<Real>, s1: &'a dyn Shape, s2: &'b dyn Shape| {
                let shape1 = s1;
                let shape2 = s2.as_shape::<S2>().ok_or(())?;
                Ok(inner(pose, shape1, shape2))
            },
        ));
    }

    pub fn add_function_known_1_and_support_map<S1: Shape + 'a>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a S1, &'b dyn SupportMap) -> bool + 'c + Send + Sync + Copy,
    ) {
        let dispatch_second = self.known_first.entry(TypeId::of::<S1>()).or_insert(vec![]);
        dispatch_second.push(Box::new(
            move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                let shape1 = s1.as_shape::<S1>().ok_or(())?;
                let shape2 = s2.as_support_map().ok_or(())?;
                Ok(inner(pose, shape1, shape2))
            },
        ));
    }

    pub fn add_function_known_2_and_support_map<S2: Shape + 'b>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a dyn SupportMap, &'b S2) -> bool + 'c + Send + Sync + Copy,
    ) {
        let dispatch_second = self.known_first.entry(TypeId::of::<S2>()).or_insert(vec![]);
        dispatch_second.push(Box::new(
            move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                let shape1 = s1.as_support_map().ok_or(())?;
                let shape2 = s2.as_shape::<S2>().ok_or(())?;
                Ok(inner(pose, shape1, shape2))
            },
        ));
    }

    pub fn add_function_all_unknown(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()>
            + 'c
            + Send
            + Sync
            + Copy,
    ) {
        self.known_none.push(Box::new(inner));
    }

    pub fn add_function_all_unknown_support_map(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a dyn SupportMap, &'b dyn SupportMap) -> bool
            + 'c
            + Send
            + Sync
            + Copy,
    ) {
        self.known_none.push(Box::new(
            move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                let shape1 = s1.as_support_map().ok_or(())?;
                let shape2 = s2.as_support_map().ok_or(())?;
                Ok(inner(pose, shape1, shape2))
            },
        ));
    }

    pub fn add_function_composite_shape_1<'d: 'c, D>(
        &mut self,
        dispatcher: &'d D,
        inner: impl Fn(&D, &Isometry<Real>, &'a dyn SimdCompositeShape, &'b dyn Shape) -> bool
            + 'c
            + Send
            + Sync
            + Copy,
    ) where
        D: ?Sized + crate::query::QueryDispatcher,
    {
        self.known_none.push(Box::new(
            move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                let shape1 = s1.as_composite_shape().ok_or(())?;
                let shape2 = s2;
                Ok(inner(dispatcher, pose, shape1, shape2))
            },
        ));
    }
    pub fn add_function_composite_shape_2<'d: 'c, D>(
        &mut self,
        dispatcher: &'d D,
        inner: impl Fn(&D, &Isometry<Real>, &'a dyn Shape, &'b dyn SimdCompositeShape) -> bool
            + 'c
            + Send
            + Sync
            + Copy,
    ) where
        D: ?Sized + crate::query::QueryDispatcher,
    {
        self.known_none.push(Box::new(
            move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                let shape1 = s1;
                let shape2 = s2.as_composite_shape().ok_or(())?;
                Ok(inner(dispatcher, pose, shape1, shape2))
            },
        ));
    }

    pub fn dispatch(
        &self,
        pose: &Isometry<Real>,
        s1: &dyn Shape,
        s2: &dyn Shape,
    ) -> Result<bool, ()> {
        let key = DispatcherTypeKey(s1.type_id(), s2.type_id());
        if let Some(func) = self.known_both.get(&key) {
            return func(pose, s1, s2);
        }

        if let Some(funcs) = self.known_first.get(&s1.type_id()) {
            for func in funcs {
                if let Ok(res) = func(pose, s1, s2) {
                    return Ok(res);
                }
            }
        }

        if let Some(funcs) = self.known_second.get(&s2.type_id()) {
            for func in funcs {
                if let Ok(res) = func(pose, s1, s2) {
                    return Ok(res);
                }
            }
        }

        for func in &self.known_none {
            if let Ok(res) = func(pose, s1, s2) {
                return Ok(res);
            }
        }

        Err(())
    }
}

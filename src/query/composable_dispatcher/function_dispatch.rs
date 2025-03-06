use core::any::TypeId;
use std::collections::HashMap;

use crate::{
    math::{Isometry, Real},
    shape::Shape,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DispatcherTypeKey(pub TypeId, pub TypeId);

pub struct FunctionDispatchFistParam<'c, 'a: 'c, 'b: 'c> {
    /// Map of functions to call when we know how to handle the first parameter.
    pub map_first_param: HashMap<TypeId, FunctionDispatchSecondParam<'c, 'a, 'b>>,
    /// Functions to attempt to call when we don't know how to handle the first parameter.
    pub dispatch_unknown_first_param: FunctionDispatchSecondParam<'c, 'a, 'b>,
}

pub struct FunctionDispatchSecondParam<'c, 'a: 'c, 'b: 'c> {
    /// Map of functions to call when we know how to handle the second parameter.
    pub map_second_param: HashMap<
        TypeId,
        Box<dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()> + 'c>,
    >,
    /// Functions to attempt to call when we don't know how to handle the second parameter.
    pub dispatch_unknown_second_param:
        Vec<Box<dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()> + 'c>>,
}

/// Transform a function that takes two concrete shapes to a function that takes two dynamic shapes.
// TODO: this would be great to have as const, or a macro.
fn to_as_shape<'c, 'a, 'b, S1: Shape + 'a, S2: Shape + 'b>(
    inner: impl Fn(&Isometry<Real>, &'a S1, &'b S2) -> bool + 'c + Copy,
) -> Box<dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()> + 'c> {
    Box::new(
        move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
            to_custom(
                |s: &'a dyn Shape| s.as_shape::<S1>(),
                |s: &'b dyn Shape| s.as_shape::<S2>(),
                inner,
            )(pose, s1, s2)
        },
    )
}
/// Transform a function that takes two concrete shapes to a function that takes two dynamic shapes.
// TODO: this would be great to have as const, or a macro.
fn to_custom<'c, 'a: 'c, 'b: 'c, ShapeIn1: 'a, ShapeIn2: 'b>(
    shape_1: impl Fn(&'a dyn Shape) -> Option<ShapeIn1> + 'a,
    shape_2: impl Fn(&'b dyn Shape) -> Option<ShapeIn2> + 'b,
    inner: impl Fn(&Isometry<Real>, ShapeIn1, ShapeIn2) -> bool + 'c,
) -> Box<dyn Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> Result<bool, ()> + 'c> {
    Box::new(
        move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
            let shape1 = shape_1(s1).ok_or(())?;
            let shape2 = shape_2(s2).ok_or(())?;
            Ok(inner(pose, shape1, shape2))
        },
    )
}

impl<'c, 'a: 'c, 'b: 'c> FunctionDispatchFistParam<'c, 'a, 'b> {
    pub fn new() -> Self {
        Self {
            map_first_param: HashMap::new(),
            dispatch_unknown_first_param: FunctionDispatchSecondParam {
                map_second_param: HashMap::new(),
                dispatch_unknown_second_param: Vec::new(),
            },
        }
    }
    pub fn add_function_known_12<S1: Shape + 'a, S2: Shape + 'b>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a S1, &'b S2) -> bool + 'c + Copy,
    ) {
        let dispatch_second = self
            .map_first_param
            .entry(TypeId::of::<S1>())
            .or_insert_with(|| FunctionDispatchSecondParam {
                map_second_param: HashMap::new(),
                dispatch_unknown_second_param: Vec::new(),
            });
        _ = dispatch_second
            .map_second_param
            .insert(TypeId::of::<S2>(), to_as_shape(inner));
    }
    pub fn add_function_known_1<S1: Shape + 'a>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a S1, &'b dyn Shape) -> bool + 'c + Copy,
    ) {
        let dispatch_second = self
            .map_first_param
            .entry(TypeId::of::<S1>())
            .or_insert_with(|| FunctionDispatchSecondParam {
                map_second_param: HashMap::new(),
                dispatch_unknown_second_param: Vec::new(),
            });
        dispatch_second.dispatch_unknown_second_param.push(Box::new(
            move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                to_custom(
                    |s: &'a dyn Shape| s.as_shape::<S1>(),
                    |s: &'b dyn Shape| Some(s),
                    inner,
                )(pose, s1, s2)
            },
        ));
    }

    pub fn add_function_known_2<S2: Shape + 'b>(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a dyn Shape, &'b S2) -> bool + 'c + Copy,
    ) {
        _ = self
            .dispatch_unknown_first_param
            .map_second_param
            .entry(TypeId::of::<S2>())
            .or_insert_with(|| {
                Box::new(
                    move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                        to_custom(
                            |s: &'a dyn Shape| Some(s),
                            |s: &'b dyn Shape| s.as_shape::<S2>(),
                            inner,
                        )(pose, s1, s2)
                    },
                )
            });
    }
    pub fn add_function_all_unknown(
        &mut self,
        inner: impl Fn(&Isometry<Real>, &'a dyn Shape, &'b dyn Shape) -> bool + 'c + Copy,
    ) {
        self.dispatch_unknown_first_param
            .dispatch_unknown_second_param
            .push(Box::new(
                move |pose: &Isometry<Real>, s1: &dyn Shape, s2: &dyn Shape| {
                    to_custom(
                        |s: &'a dyn Shape| Some(s),
                        |s: &'b dyn Shape| Some(s),
                        inner,
                    )(pose, s1, s2)
                },
            ));
    }
}

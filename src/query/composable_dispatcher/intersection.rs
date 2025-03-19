use core::any::TypeId;

use crate::{
    math::Point,
    query::{self, QueryDispatcher},
    shape::*,
};

use super::function_dispatch::FunctionDispatch;

pub fn create_intersection_dispatcher() -> FunctionDispatch {
    let mut dispatcher = FunctionDispatch::new();

    // Register intersection functions here.

    dispatcher.add_function_known_12(query::details::intersection_test_cuboid_cuboid);
    dispatcher.add_function_known_12(query::details::intersection_test_cuboid_segment);
    dispatcher.add_function_known_12(query::details::intersection_test_segment_cuboid);
    dispatcher.add_function_known_12(query::details::intersection_test_cuboid_triangle);
    dispatcher.add_function_known_12(query::details::intersection_test_triangle_cuboid);
    dispatcher.add_function_known_12(|pos12, b1, b2| {
        let p12 = Point::from(pos12.translation.vector);
        query::details::intersection_test_ball_ball(&p12, b1, b2)
    });
    dispatcher.add_function_known_1x(
        query::details::intersection_test_ball_point_query,
        vec![
            // TODO: add Tetrahedron once shape is implemented.
            TypeId::of::<Cuboid>(),
            TypeId::of::<Capsule>(),
            TypeId::of::<Triangle>(),
            TypeId::of::<Segment>(),
            TypeId::of::<Polyline>(),
            TypeId::of::<TriMesh>(),
            TypeId::of::<HeightField>(),
            #[cfg(feature = "dim2")]
            TypeId::of::<ConvexPolygon>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<ConvexPolyhedron>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<Cylinder>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<Cone>(),
            TypeId::of::<HalfSpace>(),
            // RoundShapes
            TypeId::of::<RoundShape<Cuboid>>(),
            TypeId::of::<RoundShape<Capsule>>(),
            TypeId::of::<RoundShape<Triangle>>(),
            TypeId::of::<RoundShape<Segment>>(),
            TypeId::of::<RoundShape<Polyline>>(),
            TypeId::of::<RoundShape<TriMesh>>(),
            TypeId::of::<RoundShape<HeightField>>(),
            #[cfg(feature = "dim2")]
            TypeId::of::<RoundShape<ConvexPolygon>>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<RoundShape<ConvexPolyhedron>>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<RoundShape<Cylinder>>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<RoundShape<Cone>>(),
            TypeId::of::<RoundShape<HalfSpace>>(),
        ],
    );
    dispatcher.add_function_known_1x(
        |pos12, halfspace: &HalfSpace, other| {
            query::details::intersection_test_halfspace_support_map(
                pos12,
                halfspace,
                other.as_support_map().expect("calling `as_support_map` on a non-support-map shape, make sure your type mapping is correct."),
            )
        },
        vec![
            TypeId::of::<Ball>(),
            TypeId::of::<Capsule>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<Cone>(),
            #[cfg(feature = "dim2")]
            TypeId::of::<ConvexPolygon>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<ConvexPolyhedron>(),
            TypeId::of::<Cuboid>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<Cylinder>(),
            TypeId::of::<Segment>(),
            TypeId::of::<Triangle>(),
            // RoundShapes
            TypeId::of::<RoundShape<Ball>>(),
            TypeId::of::<RoundShape<Capsule>>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<RoundShape<Cone>>(),
            #[cfg(feature = "dim2")]
            TypeId::of::<RoundShape<ConvexPolygon>>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<RoundShape<ConvexPolyhedron>>(),
            TypeId::of::<RoundShape<Cuboid>>(),
            #[cfg(feature = "dim3")]
            TypeId::of::<RoundShape<Cylinder>>(),
            TypeId::of::<RoundShape<Segment>>(),
            TypeId::of::<RoundShape<Triangle>>(),
        ],
    );
    fn intersection_test_sm_sm(
        _: &query::composable_dispatcher::ComposableQueryDispatcher,
        pos12: &crate::math::Isometry<crate::math::Real>,
        s1: &dyn Shape,
        s2: &dyn Shape,
    ) -> Result<bool, ()> {
        Ok(query::details::intersection_test_support_map_support_map(
            pos12,
            s1.as_support_map().unwrap(),
            s2.as_support_map().unwrap(),
        ))
    }
    // TODO: add support map support map
    #[cfg(feature = "dim3")]
    dispatcher.add_raw_function(
        TypeId::of::<ConvexPolyhedron>(),
        TypeId::of::<ConvexPolyhedron>(),
        |_: &query::composable_dispatcher::ComposableQueryDispatcher,
         pos12: &crate::math::Isometry<crate::math::Real>,
         s1: &dyn Shape,
         s2: &dyn Shape| {
            Ok(query::details::intersection_test_support_map_support_map(
                pos12,
                s1.as_support_map().unwrap(),
                s2.as_support_map().unwrap(),
            ))
        },
    );

    // TODO: add composite shapes
    macro_rules! shape_impls {
        ($($prefix:tt $inner:ty),*) => {
            {
                use std::collections::HashMap;
                let mut types = HashMap::new();
                $(
                    register_shape!(types, $prefix $inner);
                )*
                types
            }
        };
    }
    macro_rules! register_shape {
        ($types:expr, any $shape:ty) => {
            _ = $types.insert(TypeId::of::<$shape>(), stringify!($shape));
        };
        ($types:tt, $f:tt $shape:ty) => {
            #[cfg(feature = $f)]
            assert!($types
                .insert(TypeId::of::<$shape>(), stringify!($shape))
                .is_none());
        };
    }

    let shape_impls = shape_impls!(
        any Ball,
        any Cuboid,
        any Capsule,
        any Triangle,
        any Segment,
        any Compound,
        any Polyline,
        any TriMesh,
        any HeightField,
        "dim2" ConvexPolygon,
        "dim3" ConvexPolyhedron,
        "dim3" Cylinder,
        "dim3" Cone,
        any HalfSpace,
        // Round shapes
        any RoundShape::<Ball>,
        any RoundShape::<Cuboid>,
        any RoundShape::<Capsule>,
        any RoundShape::<Triangle>,
        any RoundShape::<Segment>,
        any RoundShape::<Compound>,
        any RoundShape::<Polyline>,
        any RoundShape::<TriMesh>,
        any RoundShape::<HeightField>,
        "dim2" RoundShape::<ConvexPolygon>,
        "dim3" RoundShape::<ConvexPolyhedron>,
        "dim3" RoundShape::<Cylinder>,
        "dim3" RoundShape::<Cone>,
        any RoundShape::<HalfSpace>
    );

    for composite_type in vec![
        TypeId::of::<Compound>(),
        TypeId::of::<Polyline>(),
        TypeId::of::<TriMesh>(),
    ] {
        for shape_type in shape_impls.iter() {
            dispatcher.add_function_dyn_dispatcher(
                |dispatcher, pos12, shape1, shape2| {
                    query::details::intersection_test_composite_shape_shape(dispatcher,
                        pos12,
                        shape1.as_composite_shape().expect("calling `as_composite_shape` on a non-support-map shape, make sure your type mapping is correct."),
                        shape2
                    )
                },
                composite_type, *shape_type.0
            );
        }
    }

    for k in dispatcher.functions.keys() {
        println!(
            "({} {})",
            shape_impls.get(&k.0).unwrap_or(&"not found"),
            shape_impls.get(&k.1).unwrap_or(&"not found")
        );
    }
    dispatcher
}

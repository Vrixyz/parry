use crate::{math::Point, query};

use super::function_dispatch::FunctionDispatchFistParam;

pub fn create_intersection_dispatcher() -> FunctionDispatchFistParam<'static, 'static, 'static> {
    let mut dispatcher = FunctionDispatchFistParam::new();

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
    dispatcher.add_function_known_12(query::details::intersection_test_segment_cuboid);
    dispatcher.add_function_known_1(query::details::intersection_test_ball_point_query);
    dispatcher.add_function_known_2(query::details::intersection_test_point_query_ball);

    /*
    // Those can't work because SupportMap is not Shape...
    // So we'd need to rely on Any...
    // but then `Any`` doesn't support casting to traits easily:
    //  - we'd need to box the shapes
    //  - or use a trait to cast to the right trait?
    dispatcher.add_function_known_1(query::details::intersection_test_halfspace_support_map);
    dispatcher.add_function_known_2(query::details::intersection_test_support_map_halfspace);

    dispatcher.add_function_all_unknown(query::details::intersection_test_support_map_support_map);
    dispatcher.add_function_known_2(query::details::intersection_test_support_map_support_map);
    */
    dispatcher
}

fn cast_params() {}

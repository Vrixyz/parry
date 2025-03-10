extern crate nalgebra as na;

#[cfg(feature = "dim2")]
mod tests2d {
    use crate::query::composable_dispatcher::intersection::create_intersection_dispatcher;
    use crate::query::composable_dispatcher::ComposableQueryDispatcher;
    use crate::query::QueryDispatcher;
    use crate::shape::{Ball, Cuboid};
    use na::{Isometry2, Vector2};

    #[test]
    fn intersection_test_dispatcher_configuration() {
        let dispatcher = create_intersection_dispatcher();
    }

    // #[test]
    // fn intersection_test() {
    //     let cuboid = Cuboid::new(Vector2::new(1.0, 1.0));
    //     let ball = Ball::new(1.0);

    //     let cuboid_pos = Isometry2::identity();
    //     let ball_pos_intersecting = Isometry2::translation(1.0, 1.0);
    //     let ball_pos_disjoint = Isometry2::translation(3.0, 3.0);

    //     let pos12 = ball_pos_intersecting.inv_mul(&cuboid_pos);
    //     assert!(ComposableQueryDispatcher
    //         .intersection_test(&pos12, &ball, &cuboid)
    //         .unwrap());

    //     let pos12 = ball_pos_disjoint.inv_mul(&cuboid_pos);
    //     assert!(!ComposableQueryDispatcher
    //         .intersection_test(&pos12, &ball, &cuboid)
    //         .unwrap());
    // }
}

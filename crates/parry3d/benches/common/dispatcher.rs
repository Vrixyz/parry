use na::Isometry3;
use parry3d::math::Isometry;
use parry3d::query::composable_dispatcher::ComposableQueryDispatcher;
use parry3d::query::{self, Unsupported};
use parry3d::query::{
    composable_dispatcher::create_intersection_dispatcher, DefaultQueryDispatcher, QueryDispatcher,
};
use parry3d::shape::{Ball, ConvexPolyhedron, Shape};
use rand::SeedableRng;
use rand_isaac::IsaacRng;
use test::Bencher;

use crate::common::{generate, unref};

#[path = "../common/macros.rs"]
#[macro_use]
mod macros;

#[bench]
fn bench_dispatcher_creation(bh: &mut Bencher) {
    let mut i = 0;
    bh.iter(|| unsafe { test::black_box(create_intersection_dispatcher()) })
}

pub fn intersection_default(
    pos12: &Isometry<f32>,
    g1: &dyn Shape,
    g2: &dyn Shape,
) -> Result<bool, Unsupported> {
    DefaultQueryDispatcher.intersection_test(pos12, g1, g2)
}

bench_free_fn!(
    bench_intersection_ball_ball,
    intersection_default,
    pos12: Isometry3<f32>,
    b1: Ball,
    b2: Ball
);

#[bench]
fn bench_intersection_ball_ball_composable(bh: &mut Bencher) {
    const LEN: usize = 1 << 7;
    let mut rng: IsaacRng = SeedableRng::seed_from_u64(0);
    let pos12: Vec<Isometry3<f32>> = (0usize..LEN).map(|_| generate(&mut rng)).collect();
    let b1: Vec<Ball> = (0usize..LEN).map(|_| generate(&mut rng)).collect();
    let b2: Vec<Ball> = (0usize..LEN).map(|_| generate(&mut rng)).collect();
    let mut i = 0;
    let dispatcher = ComposableQueryDispatcher {
        intersections: create_intersection_dispatcher(),
    };
    bh.iter(|| {
        i = (i + 1) & (LEN - 1);
        unsafe {
            test::black_box(dispatcher.intersection_test(
                unref(pos12.get_unchecked(i)),
                unref(b1.get_unchecked(i)),
                unref(b2.get_unchecked(i)),
            ))
        }
    });
}

bench_free_fn!(
    bench_intersection_convex_polyhedron_convex_polyhedron,
    intersection_default,
    pos12: Isometry3<f32>,
    b1: ConvexPolyhedron,
    b2: ConvexPolyhedron
);

#[bench]
fn bench_intersection_convex_polyhedron_convex_polyhedron_composable(bh: &mut Bencher) {
    const LEN: usize = 1 << 7;
    let mut rng: IsaacRng = SeedableRng::seed_from_u64(0);
    let pos12: Vec<Isometry3<f32>> = (0usize..LEN).map(|_| generate(&mut rng)).collect();
    let b1: Vec<ConvexPolyhedron> = (0usize..LEN).map(|_| generate(&mut rng)).collect();
    let b2: Vec<ConvexPolyhedron> = (0usize..LEN).map(|_| generate(&mut rng)).collect();
    let mut i = 0;
    let dispatcher = ComposableQueryDispatcher {
        intersections: create_intersection_dispatcher(),
    };
    bh.iter(|| {
        i = (i + 1) & (LEN - 1);
        unsafe {
            assert!(matches!(
                test::black_box(dispatcher.intersection_test(
                    unref(pos12.get_unchecked(i)),
                    unref(b1.get_unchecked(i)),
                    unref(b2.get_unchecked(i)),
                )),
                Ok(_)
            ))
        }
    });
}

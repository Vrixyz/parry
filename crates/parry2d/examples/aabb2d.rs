mod common_macroquad;

extern crate nalgebra as na;

use common_macroquad::{
    draw_point, draw_polygon, draw_polyline, lissajous_2d, lissajous_2d_with_params, mquad_from_na,
    na_from_mquad,
};
use macroquad::prelude::*;
use na::{iter, Isometry2};
use nalgebra::Point2;
use parry2d::bounding_volume::{Aabb, BoundingVolume};
use parry2d::shape::Ball;
use parry2d::transformation;
use std::f32::consts::{FRAC_PI_2, FRAC_PI_4};

const RENDER_SCALE: f32 = 30.0;

#[macroquad::main("parry2d::utils::point_in_poly2d")]
async fn main() {
    let render_pos = Vec2::new(300.0, 300.0);

    loop {
        let elapsed_time = get_time() as f32;
        let elapsed_time_slow = elapsed_time * 0.2;
        clear_background(BLACK);

        /*
         * Initialize the shapes.
         */
        let ball1 = Ball::new(0.5);
        let ball2 = Ball::new(1.0);

        let ball1_pos = na_from_mquad(lissajous_2d(elapsed_time)) * 3f32;
        let ball2_pos = Isometry2::identity();

        /*
         * Compute their axis-aligned bounding boxes.
         */
        let aabb_ball1 = ball1.aabb(&ball1_pos.into());
        let aabb_ball2 = ball2.aabb(&ball2_pos);

        // Merge the two boxes.
        let bounding_aabb = aabb_ball1.merged(&aabb_ball2);

        // Enlarge the ball2 aabb.
        let loose_aabb_ball2 = aabb_ball2.loosened(1.0);

        // Intersection and inclusion tests.
        let color = if aabb_ball1.intersects(&aabb_ball2) {
            RED
        } else {
            GREEN
        };

        assert!(bounding_aabb.contains(&aabb_ball1));
        assert!(bounding_aabb.contains(&aabb_ball2));
        assert!(loose_aabb_ball2.contains(&aabb_ball2));

        let ball1_translation =
            Vec2::new(ball1_pos.coords.x, ball1_pos.coords.y) * RENDER_SCALE + render_pos;
        draw_circle(
            ball1_translation.x,
            ball1_translation.y,
            ball1.radius * RENDER_SCALE,
            color,
        );
        let ball2_translation = Vec2::new(
            ball2_pos.translation.vector.x,
            ball2_pos.translation.vector.y,
        ) * RENDER_SCALE
            + render_pos;
        draw_circle(
            ball2_translation.x,
            ball2_translation.y,
            ball2.radius * RENDER_SCALE,
            color,
        );

        draw_aabb(aabb_ball1, render_pos, color);
        draw_aabb(aabb_ball2, render_pos, color);
        draw_aabb(bounding_aabb, render_pos, BLUE);
        draw_aabb(loose_aabb_ball2, render_pos, YELLOW);
        next_frame().await
    }
}

fn draw_aabb(aabb: Aabb, offset: Vec2, color: Color) {
    let mins = mquad_from_na(aabb.mins) * RENDER_SCALE + offset;
    let maxs = mquad_from_na(aabb.maxs) * RENDER_SCALE + offset;
    let size = maxs - mins;

    let line = vec![
        Vec2::new(mins.x, mins.y),
        Vec2::new(mins.x, maxs.y),
        Vec2::new(maxs.x, maxs.y),
        Vec2::new(maxs.x, mins.y),
        Vec2::new(mins.x, mins.y),
    ];
    let drawable_line = line
        .iter()
        .zip(line.iter().cycle().skip(1).take(line.len()))
        .map(|item| (item.0.clone(), item.1.clone()))
        .collect();
    draw_polyline(drawable_line, color);
}

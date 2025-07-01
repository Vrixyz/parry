use super::get_rr;
use core::time::Duration;
use std::boxed::Box;

use parry3d::math::Real;
use parry3d::{
    na::Vector3,
    shape::{TopologyError, TriMeshFlags},
    transformation::TriangleFacet,
};
use rerun::{LineStrip3D, LineStrips3D, Position3D, Vec3D};

/// Calls [super::init_record]
pub fn init_parry_types() {
    super::init_record(Some(Box::new(|record| {
        /// Returns [Ok] if it logged something.
        fn try_incr_time(record: &log::Record) -> Result<(), ()> {
            let key_values = record.key_values();
            let Some(time) = key_values.get(log::kv::Key::from_str("rerun_inc_time")) else {
                return Err(());
            };
            let incr = time.to_u64().unwrap_or(1) as u32;
            let time = super::TIME.fetch_add(incr, core::sync::atomic::Ordering::Relaxed);
            std::println!("{}", time + incr);
            get_rr().set_time_seconds("frame_idx", time + incr);
            Ok(())
        }
        try_incr_time(record);
        /// Returns [Ok] if it logged something.
        fn try_log_vec_triangle_facet(record: &log::Record) -> Result<(), ()> {
            let key_values = record.key_values();
            let Some(triangle) = key_values.get(log::kv::Key::from_str("rerun_vec_triangles"))
            else {
                return Err(());
            };
            let Some(points) = key_values.get(log::kv::Key::from_str("points")) else {
                return Err(());
            };
            let triangles_str = serde_json::to_string(&triangle).unwrap();
            let points_str = serde_json::to_string(&points).unwrap();
            let triangles = serde_json::from_str::<Vec<TriangleFacet>>(&triangles_str).unwrap();
            let points = serde_json::from_str::<Vec<Vector3<Real>>>(&points_str).unwrap();

            for (i, TriangleFacet { pts, .. }) in triangles.iter().enumerate() {
                let (a, b, c) = (pts[0], pts[1], pts[2]);
                let (a, b, c) = (points[a], points[b], points[c]);
                let (a, b, c) = (
                    Vec3D::new(a.x, a.y, a.z),
                    Vec3D::new(b.x, b.y, b.z),
                    Vec3D::new(c.x, c.y, c.z),
                );
                let ls = LineStrip3D::from_iter([a, b, c, a].iter());
                //triangles.push(ls);
                get_rr()
                    .log(
                        format!(
                            "triangle {i}: {}",
                            record.args().as_str().unwrap_or("triangle facets")
                        ),
                        &LineStrips3D::new([ls]),
                    )
                    .unwrap();
            }
            Ok(())
        }
        try_log_vec_triangle_facet(record);

        fn try_log_vec_points(record: &log::Record) -> Result<(), ()> {
            let key_values = record.key_values();
            let Some(points) = key_values.get(log::kv::Key::from_str("points")) else {
                return Err(());
            };
            let points_str = serde_json::to_string(&points).unwrap();

            let points = serde_json::from_str::<Vec<Vector3<Real>>>(&points_str).unwrap();
            get_rr()
                .log(
                    format!("points: {}", record.args().as_str().unwrap_or_default()),
                    &rerun::Points3D::new(
                        points
                            .iter()
                            .map(|p| Vec3D::new(p.x, p.y, p.z))
                            .collect::<Vec<_>>(),
                    ),
                )
                .unwrap();
            Ok(())
        }
        try_log_vec_points(record);
    })));
}

use super::get_rr;
use std::boxed::Box;

use parry3d::math::Real;
use parry3d::{na::Vector3, transformation::TriangleFacet};
use rerun::{LineStrip3D, LineStrips3D, Vec3D};

/// Calls [super::init_record]
pub fn init_parry_types() {
    super::init_record(Some(Box::new(|record| {
        let _ = try_incr_time(record);
        let _ = try_log_vec_triangle_facet(record);
        let _ = try_log_vec_points(record);
        let _ = try_clear(record);
        let _ = try_log_vertices_indices(record);
    })));
}

/// Returns [Ok] if it logged something.
fn try_incr_time(record: &log::Record) -> Result<(), ()> {
    let key_values = record.key_values();
    let Some(time) = key_values.get(log::kv::Key::from_str("rerun_inc_time")) else {
        return Err(());
    };
    let incr = time.to_u64().unwrap_or(1) as u32;
    let time = super::TIME.fetch_add(incr, core::sync::atomic::Ordering::Relaxed);
    std::println!("{}", time + incr);
    //get_rr().set_time_seconds("frame_idx", time + incr);
    get_rr().set_duration_secs("sim_time", time + incr);
    Ok(())
}

/// Returns [Ok] if it logged something.
fn try_log_vec_triangle_facet(record: &log::Record) -> Result<(), ()> {
    let key_values = record.key_values();
    let Some(triangle) = key_values.get(log::kv::Key::from_str("parry_vec_triangles")) else {
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
                    "{}/{i}/",
                    record.args().as_str().unwrap_or("triangle facets")
                ),
                &LineStrips3D::new([ls]),
            )
            .unwrap();
    }
    Ok(())
}

/// Returns [Ok] if it logged something.
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

/// Returns [Ok] if it logged something.
fn try_log_vertices_indices(record: &log::Record) -> Result<(), ()> {
    let key_values = record.key_values();
    let Some(vertices) = key_values.get(log::kv::Key::from_str("vertices")) else {
        return Err(());
    };
    let Some(indices) = key_values.get(log::kv::Key::from_str("indices")) else {
        return Err(());
    };
    let vertices_str = serde_json::to_string(&vertices).unwrap();
    let indices_str = serde_json::to_string(&indices).unwrap();
    let vertices = serde_json::from_str::<Vec<Vector3<Real>>>(&vertices_str).unwrap();
    let indices = serde_json::from_str::<Vec<[u32; 3]>>(&indices_str).unwrap();

    for (i, [a, b, c]) in indices.iter().enumerate() {
        let (a, b, c) = (
            vertices[*a as usize],
            vertices[*b as usize],
            vertices[*c as usize],
        );
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
                    "{}/{i}/",
                    record.args().as_str().unwrap_or("triangle facets")
                ),
                &LineStrips3D::new([ls]),
            )
            .unwrap();
    }
    Ok(())
}

/// Returns [Ok] if it logged something.
fn try_clear(record: &log::Record) -> Result<(), ()> {
    let key_values = record.key_values();
    let Some(clear_name) = key_values.get(log::kv::Key::from_str("clear")) else {
        return Err(());
    };
    let clear_name = serde_json::to_string(&clear_name).unwrap();
    let clear_name = serde_json::from_str::<String>(&clear_name).unwrap();
    println!("CLEAR: {clear_name:?}");
    get_rr()
        .log(clear_name, &rerun::Clear::recursive())
        .unwrap();
    Ok(())
}

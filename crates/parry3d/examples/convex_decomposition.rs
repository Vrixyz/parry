mod common_macroquad3d;

extern crate nalgebra as na;

use common_macroquad3d::{easy_draw_text, hue_to_rgb, mquad_mesh_from_points};
use macroquad::prelude::*;
use na::Point3;
use obj::{Obj, ObjData};
use parry3d::{
    math::Real,
    shape::{SharedShape, TriMesh, TriMeshFlags},
};

#[macroquad::main("convex_decomposition")]
async fn main() {
    /*
     * Initialize the shapes.
     */
    let Obj {
        data: ObjData {
            position, objects, ..
        },
        ..
    } = Obj::load("assets/tests/stairs.obj").unwrap();

    let bunny_mesh = TriMesh::with_flags(
        position
            .iter()
            .map(|v| Point3::new(v[0] as Real, v[1] as Real, v[2] as Real))
            .collect::<Vec<_>>(),
        objects[0].groups[0]
            .polys
            .iter()
            .map(|p| [p.0[0].0 as u32, p.0[1].0 as u32, p.0[2].0 as u32])
            .collect::<Vec<_>>(),
        TriMeshFlags::all(),
    )
    .unwrap();
    clear_background(BLACK);

    easy_draw_text("Please wait while convex decomposition is being computed...");
    next_frame().await;
    let mesh_vertices = bunny_mesh.vertices();
    let mesh_indices = bunny_mesh.indices();
    let convex_mesh = SharedShape::convex_decomposition(&mesh_vertices, &mesh_indices);
    let trimesh_convex_compound = convex_mesh.as_compound().unwrap();

    loop {
        clear_background(BLACK);
        let elapsed_time = get_time() as f32;
        let camera_pos = Vec3::new(
            19f32 * elapsed_time.sin(),
            12f32,
            19f32 * elapsed_time.cos(),
        );
        // Initialize 3D camera.
        set_camera(&Camera3D {
            position: camera_pos,
            up: Vec3::new(0f32, 1f32, 0f32),
            target: Vec3::new(0f32, 4.5f32, 0f32),
            ..Default::default()
        });

        let shapes_count = trimesh_convex_compound.shapes().len() as u32;
        for (i, s) in trimesh_convex_compound.shapes().iter().enumerate() {
            let trimesh_convex = s.1.as_convex_polyhedron().unwrap().to_trimesh();

            /*
             * Display the shapes.
             */
            let (r, g, b) = hue_to_rgb(i as f32 / shapes_count as f32);
            let mesh = mquad_mesh_from_points(
                &trimesh_convex,
                Vec3::new(1f32, 3f32, 3f32),
                Color::from_rgba(
                    (r as f32 * 255.0) as u8,
                    (g as f32 * 255.0) as u8,
                    (b as f32 * 255.0) as u8,
                    255,
                ),
            );
            draw_mesh(&mesh);
        }
        next_frame().await
    }
}

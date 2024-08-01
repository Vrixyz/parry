use na::Vector2;
use parry2d::{
    math::{Isometry, Point, Real, Vector},
    query,
    shape::{Capsule, ConvexPolygon},
};

#[test]
fn capsule_convergence() {
    let p = Point::from(Vector::y() * 5.0);

    let shape1 = Capsule::new(-p, p, 10.0);
    let mut vec = Vec::<Point<Real>>::with_capacity(3);
    vec.push(Point::<Real>::new(64.0, 507.0));
    vec.push(Point::<Real>::new(440.0, 326.0));
    vec.push(Point::<Real>::new(1072.0, 507.0));
    let shape2 = ConvexPolygon::from_convex_polyline(vec);
    let shape2 = shape2.unwrap();
    let transform1 = Isometry::new(Vector2::new(381.592, 348.491), 0.0);
    let transform2 = Isometry::new(Vector2::new(0.0, 0.0), 0.0);

    let res = query::details::contact_support_map_support_map(
        &transform1.inv_mul(&transform2),
        &shape1,
        &shape2,
        10.0,
    )
    .expect("Penetration not found.");

    if let Ok(Some(contact)) = query::contact(&transform1, &shape1, &transform2, &shape2, 1.0) {
        panic!("collision");
    } else {
        print!("no collision");
    }
}

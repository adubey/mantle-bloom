mod climate;
mod distance;
mod eustasy;
mod graph;
mod ice_step;
mod plastic;

use pyo3::prelude::*;

#[pymodule]
fn mantle_bloom_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<graph::MarchGraph>()?;
    m.add_function(wrap_pyfunction!(ice_step::step_ice, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::seeded_flood_levels, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::solve_sea_level_connected, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::below_sea_components, m)?)?;
    m.add_function(wrap_pyfunction!(eustasy::connected_ocean_mask, m)?)?;
    m.add_function(wrap_pyfunction!(climate::march_substep_prep, m)?)?;
    m.add_function(wrap_pyfunction!(climate::apply_march_abrasion, m)?)?;
    Ok(())
}

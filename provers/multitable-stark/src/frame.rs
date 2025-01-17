use crate::{table::TableView, trace::LDETraceTable};
use itertools::Itertools;
use lambdaworks_math::field::traits::{IsField, IsSubFieldOf};

/// A frame represents a collection of trace steps.
/// The collected steps are all the necessary steps for
/// all transition costraints over a trace to be evaluated.
#[derive(Clone, Debug, PartialEq)]
pub struct Frame<'t, F: IsSubFieldOf<E>, E: IsField>
where
    E: IsField,
    F: IsSubFieldOf<E>,
{
    steps: Vec<TableView<'t, F, E>>,
}

impl<'t, F: IsSubFieldOf<E>, E: IsField> Frame<'t, F, E> {
    pub fn new(steps: Vec<TableView<'t, F, E>>) -> Self {
        Self { steps }
    }

    pub fn get_evaluation_step(&self, step: usize) -> &TableView<'t, F, E> {
        &self.steps[step]
    }

    pub fn read_from_lde(
        lde_trace: &'t LDETraceTable<F, E>,
        row: usize,
        offsets: &[usize],
    ) -> Self {
        let blowup_factor = lde_trace.blowup_factor;
        let num_rows = lde_trace.num_rows();
        let step_size = lde_trace.lde_step_size;

        let lde_steps = offsets
            .iter()
            .map(|offset| {
                let initial_step_row = row + offset * step_size;
                let end_step_row = initial_step_row + step_size;
                let (table_view_main_data, table_view_aux_data) = (initial_step_row..end_step_row)
                    .step_by(blowup_factor)
                    .map(|step_row| {
                        let step_row_idx = step_row % num_rows;
                        let main_row = lde_trace.get_main_row(step_row_idx);
                        let aux_row = lde_trace.get_aux_row(step_row_idx);
                        (main_row, aux_row)
                    })
                    .unzip();

                TableView::new(table_view_main_data, table_view_aux_data)
            })
            .collect_vec();

        Frame::new(lde_steps)
    }

    pub fn read_step_from_lde(
        lde_trace: &'t LDETraceTable<F, E>,
        step: usize,
        offsets: &[usize],
    ) -> Self {
        let blowup_factor = lde_trace.blowup_factor;
        let num_rows = lde_trace.num_rows();
        let step_size = lde_trace.lde_step_size;
        let row = lde_trace.step_to_row(step);

        let lde_steps = offsets
            .iter()
            .map(|offset| {
                let initial_step_row = row + offset * step_size;
                let end_step_row = initial_step_row + step_size;
                let (table_view_main_data, table_view_aux_data) = (initial_step_row..end_step_row)
                    .step_by(blowup_factor)
                    .map(|step_row| {
                        let step_row_idx = step_row % num_rows;
                        let main_row = lde_trace.get_main_row(step_row_idx);
                        let aux_row = lde_trace.get_aux_row(step_row_idx);
                        (main_row, aux_row)
                    })
                    .unzip();

                TableView::new(table_view_main_data, table_view_aux_data)
            })
            .collect_vec();

        Frame::new(lde_steps)
    }
}

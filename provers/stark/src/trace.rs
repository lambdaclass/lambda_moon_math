use crate::table::Table;
use itertools::Itertools;
use lambdaworks_math::fft::errors::FFTError;
use lambdaworks_math::field::traits::{IsField, IsSubFieldOf};
use lambdaworks_math::{
    field::{element::FieldElement, traits::IsFFTField},
    gpu::icicle::IcicleFFT,
    polynomial::Polynomial,
    traits::ByteConversion,
};
#[cfg(feature = "parallel")]
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};

/// A two-dimensional representation of an execution trace of the STARK
/// protocol.
///
/// For the moment it is mostly a wrapper around the `Table` struct. It is a
/// layer above the raw two-dimensional table, with functionality relevant to the
/// STARK protocol, such as the step size (number of consecutive rows of the table)
/// of the computation being proven.
#[derive(Clone, Default, Debug, PartialEq, Eq)]
pub struct TraceTable<F: IsField> {
    pub table: Table<F>,
    pub step_size: usize,
    pub num_main_columns: usize,
    pub num_aux_columns: usize,
}

impl<F: IsField> TraceTable<F> {
    pub fn new(
        data: Vec<FieldElement<F>>,
        num_main_columns: usize,
        num_aux_columns: usize,
        step_size: usize,
    ) -> Self {
        let num_columns = num_main_columns + num_aux_columns;
        let table = Table::new(data, num_columns);
        Self {
            table,
            num_main_columns,
            num_aux_columns,
            step_size,
        }
    }

    pub fn from_columns(
        columns: Vec<Vec<FieldElement<F>>>,
        num_main_columns: usize,
        step_size: usize,
    ) -> Self {
        println!("COLUMNS LEN: {}", columns.len());
        println!("NUM MAIN COLUMNS: {}", num_main_columns);
        let num_aux_columns = columns.len() - num_main_columns;
        let table = Table::from_columns(columns);
        Self {
            table,
            num_main_columns,
            num_aux_columns,
            step_size,
        }
    }

    pub fn from_columns_main(columns: Vec<Vec<FieldElement<F>>>, step_size: usize) -> Self {
        let num_main_columns = columns.len();
        let num_aux_columns = 0;
        let table = Table::from_columns(columns);

        Self {
            table,
            num_main_columns,
            num_aux_columns,
            step_size,
        }
    }

    pub fn empty() -> Self {
        Self::new(Vec::new(), 0, 0, 0)
    }

    pub fn is_empty(&self) -> bool {
        self.table.width == 0
    }

    pub fn n_rows(&self) -> usize {
        self.table.height
    }

    pub fn num_steps(&self) -> usize {
        debug_assert!((self.table.height % self.step_size) == 0);
        self.table.height / self.step_size
    }

    /// Given a particular step of the computation represented on the trace,
    /// returns the row of the underlying table.
    pub fn step_to_row(&self, step: usize) -> usize {
        self.step_size * step
    }

    pub fn n_cols(&self) -> usize {
        self.table.width
    }

    pub fn rows(&self) -> Vec<Vec<FieldElement<F>>> {
        self.table.rows()
    }

    pub fn get_row(&self, row_idx: usize) -> &[FieldElement<F>] {
        self.table.get_row(row_idx)
    }

    pub fn get_row_mut(&mut self, row_idx: usize) -> &mut [FieldElement<F>] {
        self.table.get_row_mut(row_idx)
    }

    pub fn last_row(&self) -> &[FieldElement<F>] {
        self.get_row(self.n_rows() - 1)
    }

    pub fn columns(&self) -> Vec<Vec<FieldElement<F>>> {
        self.table.columns()
    }

    /// Given a slice of integer numbers representing column indexes, merge these columns into
    /// a one-dimensional vector.
    ///
    /// The particular way they are merged is not really important since this function is used to
    /// aggreagate values distributed across various columns with no importance on their ordering,
    /// such as to sort them.
    pub fn merge_columns(&self, column_indexes: &[usize]) -> Vec<FieldElement<F>> {
        let mut data = Vec::with_capacity(self.n_rows() * column_indexes.len());
        for row_index in 0..self.n_rows() {
            for column in column_indexes {
                data.push(self.table.data[row_index * self.n_cols() + column].clone());
            }
        }
        data
    }

    pub fn compute_trace_polys<S>(&self) -> Vec<Polynomial<FieldElement<F>>>
    where
        S: IsFFTField + IsSubFieldOf<F>,
        FieldElement<F>: Send + Sync + ByteConversion,
        F: IcicleFFT + IsFFTField,
    {
        let columns = self.columns();
        #[cfg(feature = "parallel")]
        let iter = columns.par_iter();
        #[cfg(not(feature = "parallel"))]
        let iter = columns.iter();

        iter.map(|col| Polynomial::interpolate_fft::<S>(col))
            .collect::<Result<Vec<Polynomial<FieldElement<F>>>, FFTError>>()
            .unwrap()
    }

    /// Given the padding length, appends the last row of the trace table
    /// that many times.
    /// This is useful for example when the desired trace length should be power
    /// of two, and only the last row is the one that can be appended without affecting
    /// the integrity of the constraints.
    pub fn pad_with_last_row(&mut self, padding_len: usize) {
        let last_row = self.last_row().to_vec();
        (0..padding_len).for_each(|_| {
            self.table.append_row(&last_row);
        })
    }

    /// Given a row index, a column index and a value, tries to set that location
    /// of the trace with the given value.
    /// The row_idx passed as argument may be greater than the max row index by 1. In this case,
    /// last row of the trace is cloned, and the value is set in that cloned row. Then, the row is
    /// appended to the end of the trace.
    pub fn set_or_extend(&mut self, row_idx: usize, col_idx: usize, value: &FieldElement<F>) {
        debug_assert!(col_idx < self.n_cols());
        // NOTE: This is not very nice, but for how the function is being used at the moment,
        // the passed `row_idx` should never be greater than `self.n_rows() + 1`. This is just
        // an integrity check for ease in the developing process, we should think a better alternative
        // in the future.
        debug_assert!(row_idx <= self.n_rows() + 1);
        if row_idx >= self.n_rows() {
            let mut last_row = self.last_row().to_vec();
            last_row[col_idx] = value.clone();
            self.table.append_row(&last_row)
        } else {
            let row = self.get_row_mut(row_idx);
            row[col_idx] = value.clone();
        }
    }
}
pub struct LDETraceTable<F, E>
where
    E: IsField,
    F: IsSubFieldOf<E>,
{
    pub(crate) main_table: Table<F>,
    pub(crate) aux_table: Table<E>,
    pub(crate) lde_step_size: usize,
    pub(crate) blowup_factor: usize,
}

impl<F, E> LDETraceTable<F, E>
where
    E: IsField,
    F: IsSubFieldOf<E>,
{
    pub fn new(
        main_data: Vec<FieldElement<F>>,
        aux_data: Vec<FieldElement<E>>,
        n_columns: usize,
        trace_step_size: usize,
        blowup_factor: usize,
    ) -> Self {
        let main_table = Table::new(main_data, n_columns);
        let aux_table = Table::new(aux_data, n_columns);
        let lde_step_size = trace_step_size * blowup_factor;

        Self {
            main_table,
            aux_table,
            lde_step_size,
            blowup_factor,
        }
    }

    pub fn from_columns(
        main_columns: Vec<Vec<FieldElement<F>>>,
        aux_columns: Vec<Vec<FieldElement<E>>>,
        trace_step_size: usize,
        blowup_factor: usize,
    ) -> Self {
        let main_table = Table::from_columns(main_columns);
        let aux_table = Table::from_columns(aux_columns);
        let lde_step_size = trace_step_size * blowup_factor;

        Self {
            main_table,
            aux_table,
            lde_step_size,
            blowup_factor,
        }
    }

    pub fn num_cols(&self) -> usize {
        self.main_table.width + self.aux_table.width
    }

    pub fn num_rows(&self) -> usize {
        self.main_table.height
    }

    pub fn get_main_row(&self, row_idx: usize) -> &[FieldElement<F>] {
        self.main_table.get_row(row_idx)
    }

    pub fn get_aux_row(&self, row_idx: usize) -> &[FieldElement<E>] {
        self.aux_table.get_row(row_idx)
    }

    pub fn get_main(&self, row: usize, col: usize) -> &FieldElement<F> {
        self.main_table.get(row, col)
    }

    pub fn get_aux(&self, row: usize, col: usize) -> &FieldElement<E> {
        self.aux_table.get(row, col)
    }

    pub fn num_steps(&self) -> usize {
        debug_assert!((self.main_table.height % self.lde_step_size) == 0);
        self.main_table.height / self.lde_step_size
    }

    pub fn step_to_row(&self, step: usize) -> usize {
        self.lde_step_size * step
    }
}

/// Given a slice of trace polynomials, an evaluation point `x`, the frame offsets
/// corresponding to the computation of the transitions, and a primitive root,
/// outputs the trace evaluations of each trace polynomial over the values used to
/// compute a transition.
/// Example: For a simple Fibonacci computation, if t(x) is the trace polynomial of
/// the computation, this will output evaluations t(x), t(g * x), t(g^2 * z).
pub fn get_trace_evaluations<F, E>(
    main_trace_polys: &[Polynomial<FieldElement<F>>],
    aux_trace_polys: &[Polynomial<FieldElement<E>>],
    x: &FieldElement<E>,
    frame_offsets: &[usize],
    primitive_root: &FieldElement<F>,
    step_size: usize,
) -> Table<E>
where
    F: IsSubFieldOf<E>,
    E: IsField,
{
    let evaluation_points = frame_offsets
        .iter()
        .flat_map(|offset| {
            let exponents_range_start = offset * step_size;
            let exponents_range_end = (offset + 1) * step_size;
            (exponents_range_start..exponents_range_end).collect_vec()
        })
        .map(|exponent| primitive_root.pow(exponent) * x)
        .collect_vec();

    let main_evaluations = evaluation_points
        .iter()
        .map(|eval_point| {
            main_trace_polys
                .iter()
                .map(|main_poly| main_poly.evaluate(eval_point))
                .collect_vec()
        })
        .collect_vec();

    let aux_evaluations = evaluation_points
        .iter()
        .map(|eval_point| {
            aux_trace_polys
                .iter()
                .map(|aux_poly| aux_poly.evaluate(eval_point))
                .collect_vec()
        })
        .collect_vec();

    debug_assert_eq!(main_evaluations.len(), aux_evaluations.len());
    let mut main_evaluations = main_evaluations;
    let mut table_data = Vec::new();
    for (main_row, aux_row) in main_evaluations.iter_mut().zip(aux_evaluations) {
        main_row.extend_from_slice(&aux_row);
        table_data.extend_from_slice(main_row);
    }

    let main_trace_width = main_trace_polys.len();
    let aux_trace_width = aux_trace_polys.len();
    let table_width = main_trace_width + aux_trace_width;

    Table::new(table_data, table_width)
}

pub fn columns2rows<F: IsField>(columns: Vec<Vec<FieldElement<F>>>) -> Vec<Vec<FieldElement<F>>> {
    let num_rows = columns[0].len();
    let num_cols = columns.len();

    (0..num_rows)
        .map(|row_index| {
            (0..num_cols)
                .map(|col_index| columns[col_index][row_index].clone())
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::TraceTable;
    use lambdaworks_math::field::{element::FieldElement, fields::u64_prime_field::F17};
    type FE = FieldElement<F17>;

    #[test]
    fn test_cols() {
        let col_1 = vec![FE::from(1), FE::from(2), FE::from(5), FE::from(13)];
        let col_2 = vec![FE::from(1), FE::from(3), FE::from(8), FE::from(21)];

        let trace_table = TraceTable::from_columns(vec![col_1.clone(), col_2.clone()], 2, 1);
        let res_cols = trace_table.columns();

        assert_eq!(res_cols, vec![col_1, col_2]);
    }
}

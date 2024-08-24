use criterion::measurement::{Measurement, ValueFormatter};
use criterion::Throughput;
use humansize::{format_size, DECIMAL};

pub struct ProofSize(usize);

impl Measurement for ProofSize {
    type Intermediate = ();
    type Value = usize;

    fn start(&self) -> Self::Intermediate {
        ()
    }

    fn end(&self, _i: Self::Intermediate) -> Self::Value {
        self.0
    }

    fn add(&self, v1: &Self::Value, v2: &Self::Value) -> Self::Value {
        v1 + v2
    }

    fn zero(&self) -> Self::Value {
        0
    }

    fn to_f64(&self, value: &Self::Value) -> f64 {
        *value as f64
    }

    fn formatter(&self) -> &dyn ValueFormatter {
        &ProofSizeFormatter {}
    }
}

pub struct ProofSizeFormatter;
impl ValueFormatter for ProofSizeFormatter {
    fn format_value(&self, value: f64) -> String {
        format_size(value as u64, DECIMAL)
    }

    fn scale_values(&self, typical_value: f64, values: &mut [f64]) -> &'static str {
        // Stupid hack to allow us to return a &'static str
        match format_size(typical_value as u64, DECIMAL)
            .split(' ')
            .last()
            .unwrap()
            .to_owned()
            .as_str()
        {
            "B" => "bytes",
            "KB" => "KB",
            "MB" => "MB",
            "GB" => "GB",
            "TB" => "TB",
            &_ => unreachable!(),
        }
    }

    fn scale_throughputs(
        &self,
        typical_value: f64,
        throughput: &Throughput,
        values: &mut [f64],
    ) -> &'static str {
        // Stupid hack to allow us to return a &'static str
        match format_size(typical_value as u64, DECIMAL)
            .split(' ')
            .last()
            .unwrap()
            .to_owned()
            .as_str()
        {
            "B" => "bytes",
            "KB" => "KB",
            "MB" => "MB",
            "GB" => "GB",
            "TB" => "TB",
            &_ => unreachable!(),
        }
    }

    fn scale_for_machines(&self, values: &mut [f64]) -> &'static str {
        "bytes"
    }
}

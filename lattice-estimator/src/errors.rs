use std::fmt;
use std::fmt::{Debug, Display};

pub struct LatticeEstimatorError {
    pub(crate) message: String,
}

impl Display for LatticeEstimatorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "LatticeEstimatorError: {}", self.message)
    }
}

impl Debug for LatticeEstimatorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "LatticeEstimatorError: {}", self.message)
    }
}

impl From<String> for LatticeEstimatorError {
    fn from(message: String) -> Self {
        LatticeEstimatorError { message }
    }
}


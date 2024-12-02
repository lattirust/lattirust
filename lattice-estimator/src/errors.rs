use std::fmt;
use std::fmt::{Debug, Display};

#[derive(Debug)]
pub enum LatticeEstimatorError {
    /// Invalid parameter error with details
    InvalidParameter {
        param_name: String,
        reason: String,
    },
    /// Computation-related error
    ComputationError {
        message: String,
    },
    /// Configuration error with a specific message
    ConfigurationError {
        message: String,
    },
}

impl Display for LatticeEstimatorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            LatticeEstimatorError::InvalidParameter { param_name, reason } => {
                write!(f, "Invalid Parameter '{}': {}", param_name, reason)
            }
            LatticeEstimatorError::ComputationError { message } => {
                write!(f, "Computation Error: {}", message)
            }
            LatticeEstimatorError::ConfigurationError { message } => {
                write!(f, "Configuration Error: {}", message)
            }
        }
    }
}

impl std::error::Error for LatticeEstimatorError {}

impl From<(&str, &str)> for LatticeEstimatorError {
    fn from((param_name, reason): (&str, &str)) -> Self {
        LatticeEstimatorError::InvalidParameter {
            param_name: param_name.to_string(),
            reason: reason.to_string(),
        }
    }
}

impl From<&str> for LatticeEstimatorError {
    fn from(message: &str) -> Self {
        LatticeEstimatorError::ComputationError {
            message: message.to_string(),
        }
    }
}

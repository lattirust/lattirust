use std::fmt;
use std::fmt::Display;

#[derive(Clone, Copy)]
pub enum Norm {
    L2,
    Linf,
}

impl Display for Norm {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Norm::L2 => write!(f, "L2"),
            Norm::Linf => write!(f, "Linf"),
        }
    }
}
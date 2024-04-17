#[macro_export]
macro_rules! check {
    ($cond: expr) => {
        if !($cond) {
            return Err(ProofError::InvalidProof);
        }
    };
    ($cond: expr , $arg: expr) => {
        debug!("check!({:?}, {}) => {}", $cond, $arg, $cond);
        if !($cond) {
            return Err(ProofError::InvalidProof);
        }
    };
}

#[macro_export]
macro_rules! check_eq {
    ($a: expr, $b : expr) => {
        debug!("check_eq!(\n\t{:?}, \n\t{:?}\n) => {}", $a, $b, $a == $b);
        check!($a == $b);
    };
    ($a: expr, $b: expr, $arg: expr) => {
        debug!("check_eq!(\n\t{:?}, \n\t{:?}, \n\t{}\n) => {}", $a, $b, $arg, $a == $b);
        check!($a == $b);
    };
}

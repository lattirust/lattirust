#[macro_export]
macro_rules! check {
    ($ cond: expr) => {
        {
            if !($cond) {
                return Err(ProofError::InvalidProof);
            }
        }
    };
    ( $ cond : expr , $ ( $ arg : tt ) + ) => {
        {
            if !($cond) {
                return Err(ProofError::InvalidProof);
            }
        }
    };
}

#[macro_export]
macro_rules! check_eq {
    ( $ a : expr , $ b : expr ) => { check!($a == $b) };
    ( $ a : expr , $ b : expr , $ ( $ arg : tt ) + ) => { check!($a == $b, $($arg)+) };
}
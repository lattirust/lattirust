#[cfg(test)]
mod test {
    use lattirust_arithmetic::ntt::ntt_modulus;
    use crate::bfv::{rand_tuple, Plaintext, SecretKey};

    const N: usize = 128;
    const P: u64 = ntt_modulus::<N>(15);
    const Q: u64 = P * 5;

    #[test]
    fn zero_message() {
        let sk: SecretKey<Q, P, N> = SecretKey::new();
        let pk = sk.pk_gen();
        let ptxt = Plaintext::zero();
        let ctxt = pk.encrypt(&ptxt, rand_tuple::<Q, N>(None));
        
        let act = sk.decrypt(ctxt).poly;
        let exp = ptxt.poly;
        assert_eq!(act, exp);
    }

    #[test]
    fn rand_message() {
        let sk: SecretKey<Q, P, N> = SecretKey::new();
        let pk = sk.pk_gen();
        let ptxt = Plaintext::rand_message();
        let ctxt = pk.encrypt(&ptxt, rand_tuple::<Q, N>(None));
        
        let act = sk.decrypt(ctxt).poly;
        let exp = ptxt.poly;
        assert_eq!(act, exp);
    }
}
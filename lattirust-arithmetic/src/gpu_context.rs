#![allow(non_snake_case)]
use ark_ff::PrimeField;
use icicle_core::ntt::{self, initialize_domain};
use icicle_core::traits::FieldImpl;
use icicle_runtime::memory::{DeviceVec, HostSlice};
use icicle_runtime::{eIcicleError, runtime, Device};
use icicle_babybear::field::ScalarField as BabybearField;

use crate::ring::Zq;


use log::{info, warn, debug, error};
use pretty_env_logger::env_logger;
use std::env;


fn init() {
    let _ = env_logger::builder().is_test(false).try_init();
}

#[derive(Debug)]
pub enum Error {
    ErrNone,
    IcicleBackendNotFound,
    GpuNotAvailable,
    CopyFailed,
    InvalidInput
}

impl Error {
    /// Converts the enum to a `u8` representation
    pub fn val(&self) -> u8 {
        match self {
            Error::ErrNone => 0,
            Error::IcicleBackendNotFound => 1,
            Error::GpuNotAvailable => 2,
            Error::CopyFailed => 3,
            Error::InvalidInput => 4,
        }
    }
}

use Error::*;


// -----------------------------------------------------------------------------------------------
// GPU State

const ID : i32 = 0;
const DEVICE_NAME : &str = "CUDA";

// -----------------------------------------------------------------------------------------------
// Utils

fn load_backend(key: &str) -> Result<(), eIcicleError> {
    match env::var(key) {
        Ok(_value) => {
            return runtime::load_backend(&_value);
        },
        Err(_e) => {
            return Err(eIcicleError::UnknownError);
        },
    }
}


// Load backend and set device
pub fn try_load_and_set_GPU_backend_device() -> Error {
    init();

    let _res = load_backend("ICICLE_BACKEND_INSTALL_DIR");
    if _res.is_err() {
        warn!("Failed to set device to GPU, defaulting to CPU implementation");
        return IcicleBackendNotFound;
    }
    
    // Register GPU
    let device = Device::new(DEVICE_NAME, ID); 
    let res : Result<(), eIcicleError> = icicle_runtime::set_device(&device);
    if res.is_err() {
        warn!("Failed to set device to GPU, defaulting to CPU implementation");
        return GpuNotAvailable;
    }
    return ErrNone;
}

/// Converts a slice of Zq elements to a Vec of BabybearField elements
pub fn convert_to_babybear<const Q: u64, const N: usize>(
    input: [Zq<Q>; N],
) -> Vec<BabybearField> {
    input
        .iter()
        .map(|zq| BabybearField::from([zq.into_bigint().0[0] as u32])) 
        .collect()
}

/// Converts a slice of BabybearField elements to a Vec of Zq elements
pub fn convert_to_Zq<const Q: u64, const N: usize>(
    input: Vec<BabybearField>,
) -> Result<[Zq<Q>; N], Error> {
    if input.len() != N {
        error!("Input length does not match expected output length");
        return Err(InvalidInput);
    }

    let mut output = [Zq::<Q>::default(); N];
    for (i, babybear) in input.iter().enumerate() {
        let bytes = babybear.to_bytes_le();

        let raw_value = bytes
            .iter()
            .enumerate()
            .fold(0u64, |acc, (index, &byte)| acc | ((byte as u64) << (index * 8)));

        // Create Zq<Q> with the extracted value
        output[i] = Zq::<Q>::new(raw_value.into());
    }
    Ok(output)
}

/// Copies data from host to device
pub fn copy_from_host(input: &mut DeviceVec<BabybearField>, data: Vec<BabybearField>) -> Result<(), Error> {

    let res = input.copy_from_host(HostSlice::from_slice(&data[..]));
    if res.is_err(){
        error!("Failed to copy data from host to device");
        return Err(CopyFailed);
    }
    Ok(())
}


/// Copies data from device to host
pub fn copy_to_host(host_results: &mut Vec<BabybearField>, output: DeviceVec<BabybearField>) -> Result<(), Error> {

    let res = output.copy_to_host(HostSlice::from_mut_slice(&mut host_results[..]));
    if res.is_err(){
        error!("Failed to copy data from device to host");
        return Err(CopyFailed);
    }
    Ok(())
}

// #######################################################################################################################################
// NTT

/// Allocates memory on the device and initializes the NTT context
pub fn init_ntt_context_on_device(input_size: usize) -> Result<(DeviceVec<BabybearField>, DeviceVec<BabybearField>), eIcicleError> {
    info!("Initializing NTT context on device: Memory allocation and NTT domain initialization");

    let ntt_input = DeviceVec::<BabybearField>::device_malloc(input_size);
    let ntt_results = DeviceVec::<BabybearField>::device_malloc(input_size);

    if let (Err(_), _) | (_, Err(_)) = (&ntt_input, &ntt_results) {
        error!("Failed to allocate memory on device");
        return Err(eIcicleError::OutOfMemory);
    }

    // Initialize the NTT domain
    if let Err(_) = initialize_ntt_domain(input_size) {
        error!("Failed to initialize NTT domain");
        return Err(eIcicleError::UnknownError);
    }

    Ok((ntt_input.unwrap(), ntt_results.unwrap()))
}

pub fn get_default_ntt_config() -> ntt::NTTConfig<BabybearField> {
    ntt::NTTConfig::<BabybearField>::default()
}


fn initialize_ntt_domain(t:usize) -> Result<(), eIcicleError> {

    let res = t.try_into();
    if res.is_err() {
        return Err(eIcicleError::UnknownError);
    }
    let size = res.unwrap();
    initialize_domain(
        ntt::get_root_of_unity::<BabybearField>(size),
        &ntt::NTTInitDomainConfig::default(),
    )
}

// #######################################################################################################################################

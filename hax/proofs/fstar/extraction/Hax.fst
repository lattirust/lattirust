module Hax
#set-options "--fuel 0 --ifuel 1 --z3rlimit 15"
open Core
open FStar.Mul

let _ =
  (* This module has implicit dependencies, here we make them explicit. *)
  (* The implicit dependencies arise from typeclasses instances. *)
  let open Num_traits.Identities in
  let open Num_traits.Pow in
  let open Rayon.Iter in
  let open Rayon.Iter.Enumerate in
  let open Rayon.Iter.Map in
  let open Rayon.Slice in
  let open Rayon.Vec in
  let open Rounded_div in
  ()

/// Given a vector of vectors `v`, pads each row to the same length $l = \max_i \texttt{v}\[i\]\texttt{.len()}$ and transposes the result. The output is a Vec of Vec of dimensionts `l` times `v.len()`.
let recompose
      (#v_A #v_B: Type0)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i2: Core.Ops.Arith.t_Mul v_A v_B)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i3: Core.Marker.t_Copy v_A)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i4: Core.Iter.Traits.Accum.t_Sum v_A v_A)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i5: Core.Marker.t_Send v_A)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i6: Core.Marker.t_Sync v_A)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i7: Core.Clone.t_Clone v_B)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i8: Core.Ops.Arith.t_Mul v_B v_B)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i9: Core.Ops.Arith.t_Add v_B v_B)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i10: Num_traits.Pow.t_Pow v_B (t_Array u64 (sz 1)))
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i11: Core.Marker.t_Send v_B)
      (#[FStar.Tactics.Typeclasses.tcresolve ()] i12: Core.Marker.t_Sync v_B)
      (v: Alloc.Vec.t_Vec v_A Alloc.Alloc.t_Global)
      (b: v_B)
    : v_A =
  Rayon.Iter.f_sum #(Rayon.Iter.Map.t_Map
        (Rayon.Iter.Enumerate.t_Enumerate (Rayon.Slice.t_Iter v_A)) ((usize & v_A) -> v_A))
    #v_A
    (Rayon.Iter.f_map #(Rayon.Iter.Enumerate.t_Enumerate (Rayon.Slice.t_Iter v_A))
        #v_A
        (Rayon.Iter.f_enumerate #(Rayon.Slice.t_Iter v_A)
            (Rayon.Iter.f_par_iter #(Alloc.Vec.t_Vec v_A Alloc.Alloc.t_Global) v
              <:
              Rayon.Slice.t_Iter v_A)
          <:
          Rayon.Iter.Enumerate.t_Enumerate (Rayon.Slice.t_Iter v_A))
        (fun temp_0_ ->
            let i, vv_i:(usize & v_A) = temp_0_ in
            vv_i *!
            (Num_traits.Pow.f_pow #v_B
                #(t_Array u64 (sz 1))
                (Core.Clone.f_clone #v_B b <: v_B)
                (let list = [cast (i <: usize) <: u64] in
                  FStar.Pervasives.assert_norm (Prims.eq2 (List.Tot.length list) 1);
                  Rust_primitives.Hax.array_of_list 1 list)
              <:
              v_B)
            <:
            v_A)
      <:
      Rayon.Iter.Map.t_Map (Rayon.Iter.Enumerate.t_Enumerate (Rayon.Slice.t_Iter v_A))
        ((usize & v_A) -> v_A))

/// Returns the floor of the logarithm of `x` in base `base`.
let floor_log (x base: u128) : u128 =
  let _:Prims.unit =
    if ~.(base >. pub_u128 1 <: bool)
    then
      Rust_primitives.Hax.never_to_any (Core.Panicking.panic "assertion failed: base > 1"
          <:
          Rust_primitives.Hax.t_Never)
  in
  let _:Prims.unit =
    if ~.(x >=. pub_u128 1 <: bool)
    then
      Rust_primitives.Hax.never_to_any (Core.Panicking.panic "assertion failed: x >= 1"
          <:
          Rust_primitives.Hax.t_Never)
  in
  let floor_log:u128 = pub_u128 0 in
  let floor_log, x:(u128 & u128) =
    Core.Iter.Traits.Iterator.f_fold (Core.Iter.Traits.Collect.f_into_iter #(Core.Ops.Range.t_RangeFrom
            i32)
          ({ Core.Ops.Range.f_start = 0l } <: Core.Ops.Range.t_RangeFrom i32)
        <:
        Core.Ops.Range.t_RangeFrom i32)
      (floor_log, x <: (u128 & u128))
      (fun temp_0_ temp_1_ ->
          let floor_log, x:(u128 & u128) = temp_0_ in
          let _:i32 = temp_1_ in
          if x >=. base <: bool
          then
            let x:u128 = x /! base in
            let floor_log:u128 = floor_log +! pub_u128 1 in
            floor_log, x <: (u128 & u128)
          else floor_log, x <: (u128 & u128))
  in
  floor_log

/// Returns the maximum number of terms in the balanced decomposition in basis `b` of any `x` with $|\textt{x}| \leq \textt{max}$.
let balanced_decomposition_max_length (b max: u128) : usize =
  if max =. pub_u128 0
  then sz 0
  else ((cast (floor_log max b <: u128) <: usize) +! sz 1 <: usize) +! sz 1

/// Returns the balanced decomposition of a slice as a Vec of Vecs.
/// # Arguments
/// * `v`: input element
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
/// # Output
/// Returns `d`, the decomposition in basis `b` as a Vec of size `decomp_size`, i.e.,
/// $\texttt{v}\[i\] = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[j\]$ and $|\texttt{d}\[j\]| \leq \left\lfloor\frac{\texttt{b}}{2}\right\rfloor$.
let decompose_balanced
      (v: i128)
      (b: u128)
      (padding_size: Core.Option.t_Option usize)
      (length: usize)
    : Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global =
  let _:Prims.unit =
    if
      ~.((~.(Num_traits.Identities.f_is_zero #u128 b <: bool) <: bool) &&
        (~.(Num_traits.Identities.f_is_one #u128 b <: bool) <: bool))
    then
      Rust_primitives.Hax.never_to_any (Core.Panicking.panic_fmt (Core.Fmt.impl_2__new_const (Rust_primitives.unsize
                    (let list = ["cannot decompose in basis 0 or 1"] in
                      FStar.Pervasives.assert_norm (Prims.eq2 (List.Tot.length list) 1);
                      Rust_primitives.Hax.array_of_list 1 list)
                  <:
                  t_Slice string)
              <:
              Core.Fmt.t_Arguments)
          <:
          Rust_primitives.Hax.t_Never)
  in
  let _:Prims.unit =
    match b %! pub_u128 2, pub_u128 0 <: (u128 & u128) with
    | left_val, right_val ->
      if ~.(left_val =. right_val <: bool)
      then
        let kind:Core.Panicking.t_AssertKind =
          Core.Panicking.AssertKind_Eq <: Core.Panicking.t_AssertKind
        in
        Rust_primitives.Hax.never_to_any (Core.Panicking.assert_failed #u128
              #u128
              kind
              left_val
              right_val
              (Core.Option.Option_Some
                (Core.Fmt.impl_2__new_const (Rust_primitives.unsize (let list =
                            ["decomposition basis must be even"]
                          in
                          FStar.Pervasives.assert_norm (Prims.eq2 (List.Tot.length list) 1);
                          Rust_primitives.Hax.array_of_list 1 list)
                      <:
                      t_Slice string)
                  <:
                  Core.Fmt.t_Arguments)
                <:
                Core.Option.t_Option Core.Fmt.t_Arguments)
            <:
            Rust_primitives.Hax.t_Never)
  in
  let b_half_floor:u128 = Core.Num.impl__u128__div_euclid b (pub_u128 2) in
  let b:i128 = cast (b <: u128) <: i128 in
  let decomp_bal_signed:Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global = Alloc.Vec.impl__new #i128 () in
  let curr:i128 = v in
  let curr, decomp_bal_signed:(i128 & Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global) =
    Core.Iter.Traits.Iterator.f_fold (Core.Iter.Traits.Collect.f_into_iter #(Core.Ops.Range.t_Range
            usize)
          ({ Core.Ops.Range.f_start = sz 0; Core.Ops.Range.f_end = length }
            <:
            Core.Ops.Range.t_Range usize)
        <:
        Core.Ops.Range.t_Range usize)
      (curr, decomp_bal_signed <: (i128 & Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global))
      (fun temp_0_ temp_1_ ->
          let curr, decomp_bal_signed:(i128 & Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global) =
            temp_0_
          in
          let _:usize = temp_1_ in
          let rem:i128 = curr %! b in
          if (cast (Core.Num.impl__i128__abs rem <: i128) <: u128) <=. b_half_floor
          then
            let decomp_bal_signed:Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global =
              Alloc.Vec.impl_1__push #i128 #Alloc.Alloc.t_Global decomp_bal_signed rem
            in
            let curr:i128 = curr /! b in
            curr, decomp_bal_signed <: (i128 & Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global)
          else
            let decomp_bal_signed:Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global =
              if rem <. pub_i128 0
              then
                let decomp_bal_signed:Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global =
                  Alloc.Vec.impl_1__push #i128
                    #Alloc.Alloc.t_Global
                    decomp_bal_signed
                    (rem +! b <: i128)
                in
                decomp_bal_signed
              else
                let decomp_bal_signed:Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global =
                  Alloc.Vec.impl_1__push #i128
                    #Alloc.Alloc.t_Global
                    decomp_bal_signed
                    (rem -! b <: i128)
                in
                decomp_bal_signed
            in
            let carry:i128 = Rounded_div.f_rounded_div #i128 #i128 rem b in
            let curr:i128 = (curr /! b <: i128) +! carry in
            curr, decomp_bal_signed <: (i128 & Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global))
  in
  let decomp_bal_signed:Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global =
    match padding_size with
    | Core.Option.Option_Some padding_size ->
      let _:Prims.unit =
        if
          ~.((Alloc.Vec.impl_1__len #i128 #Alloc.Alloc.t_Global decomp_bal_signed <: usize) <=.
            padding_size
            <:
            bool)
        then
          Rust_primitives.Hax.never_to_any (Core.Panicking.panic_fmt (Core.Fmt.impl_2__new_v1 (Rust_primitives.unsize
                        (let list = ["padding_size = "; " must be at least decomp_bal.len() = "] in
                          FStar.Pervasives.assert_norm (Prims.eq2 (List.Tot.length list) 2);
                          Rust_primitives.Hax.array_of_list 2 list)
                      <:
                      t_Slice string)
                    (Rust_primitives.unsize (let list =
                            [
                              Core.Fmt.Rt.impl_1__new_display #usize padding_size
                              <:
                              Core.Fmt.Rt.t_Argument;
                              Core.Fmt.Rt.impl_1__new_display #usize
                                (Alloc.Vec.impl_1__len #i128 #Alloc.Alloc.t_Global decomp_bal_signed
                                  <:
                                  usize)
                              <:
                              Core.Fmt.Rt.t_Argument
                            ]
                          in
                          FStar.Pervasives.assert_norm (Prims.eq2 (List.Tot.length list) 2);
                          Rust_primitives.Hax.array_of_list 2 list)
                      <:
                      t_Slice Core.Fmt.Rt.t_Argument)
                  <:
                  Core.Fmt.t_Arguments)
              <:
              Rust_primitives.Hax.t_Never)
      in
      let decomp_bal_signed:Alloc.Vec.t_Vec i128 Alloc.Alloc.t_Global =
        Alloc.Vec.impl_2__resize #i128
          #Alloc.Alloc.t_Global
          decomp_bal_signed
          padding_size
          (pub_i128 0)
      in
      decomp_bal_signed
    | _ -> decomp_bal_signed
  in
  decomp_bal_signed

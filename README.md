# FIPS-203

This is an implementation of the FIPS-203 post-quantum key encapsulation algorithm in pure rust.

> [!WARNING]
> DANGER WILL ROBINSON! Very much not yet complete. Use at your own peril ..

Spec: https://csrc.nist.gov/pubs/fips/203/final \
Static test vectors: https://github.com/usnistgov/ACVP-Server/tree/master/gen-val/json-files \
Static test vectors with intermediate values: https://github.com/C2SP/CCTV \
Reference implementation for FIPS203/KEM: https://github.com/pq-crystals/kyber.git \

## Running long-running tests

There are some very long-running tests (in particular, the reference test vectors for ML-KEM take some minutes to run through).

As a result, these tests are marked `#[ignore]` . To run them:

```
git clone git@github.com:rrw1000/fips-vectors
cd fips-vectors/kem
bunzip2 *.bz2
cd ../..
cargo test vectors -- --ignored
```

I did try `#cfg()` , but this requires `profile-rustflags` support in cargo, which is not yet (1.86) stable.

## A note on integer types and array lengths

We mostly use `i32` throughout (at least, throughout `dsa/`). This is
cheating, because sometimes we do put -ve numbers in vectors and in
other places, and this does mean there is the possibility that we will
pass -ve numbers to functions like `integer_to_bits()` which cannot
handle them.

What we should do is use `u32` where we want only +ve numbers and
`i32` everywhere else, but this would mean a lot of tediously
translating vectors from `[i32;256]` to `[u32;256]`. So we don't.

A similar argument applies to `&[u8]` vs `&[u8;32]`.

Perhaps we should be more strict .. either way, beware.

## Things to do

Things to do in future:

  * Polynomials should have their own type aliases.
  * Profile and speed up.
  * Accumulate currently doesn't, really - it copies; we should make it reuse storage if faster (it may not - in fact, probably isn't - faster)
  * BytesToBits and BitsToBytes could probably be faster (but does this make a difference in practice?)
  * We need to make as much of this as possible constant time.
  

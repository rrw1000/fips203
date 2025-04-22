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

## Things to do

Things to do in future:

  * Polynomials should have their own type aliases.
  * Profile and speed up.
  * Accumulate currently doesn't, really - it copies; we should make it reuse storage if faster (it may not - in fact, probably isn't - faster)
  * BytesToBits and BitsToBytes could probably be faster (but does this make a difference in practice?)
  K

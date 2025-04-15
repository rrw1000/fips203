# FIPS-203

This is an implementation of the FIPS-203 post-quantum key encapsulation algorithm in pure rust.

> [!WARNING]
> DANGER WILL ROBINSON! Very much not yet complete. Use at your own peril ..

Spec: https://csrc.nist.gov/pubs/fips/203/final
Static test vectors: https://github.com/usnistgov/ACVP-Server/tree/master/gen-val/json-files
Static test vectors with intermediate values: https://github.com/C2SP/CCTV

## Running long-running tests

There are some very long-running tests (in particular, the reference test vectors for ML-KEM take some minutes to run through).

As a result, these tests are marked `#[ignore]` . To run them:

```
$ cargo test vectors -- --ignored
```

I did try `#cfg()` , but this requires `profile-rustflags` support in cargo, which is not yet (1.86) stable.


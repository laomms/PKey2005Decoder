#pragma once
#include "pch.h"
#include <vector>
#include <string>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <cmath>
#include "BigInteger.h"
#include "Montgomery.h"
#include "TatePairing.h"
#include "H1Search.h"

// -------------------------------------------------------------------------
// Public Key Definition
// -------------------------------------------------------------------------
struct PubKey {
    int SizeModulus = 0;
    int SizeOrder = 0;
    int ExtDeg1 = 0;
    int ExtDeg2 = 0;
    std::vector<unsigned char> H1Bases;
    BigInteger Modulus = 0;
    BigInteger Order = 0;
    std::vector<unsigned char> K3Minpoly;
    std::vector<unsigned char> K6Minpoly;
    BigInteger CurveA = 0;
    BigInteger CurveB = 0;
    std::vector<std::vector<Fp3>> Points; // Each entry contains [x, y]
    std::vector<Fp3> PairingVal;          // K6 = [c0, c1] represented as Fp3 array/parts
};

// -------------------------------------------------------------------------
// Public Key Parser
// -------------------------------------------------------------------------
class PubKeyParser {
public:
    static PubKey Parse(const std::vector<unsigned char>& data) {
        if (data.size() < 59) {
            throw std::runtime_error("Public key data too short.");
        }

        auto U32 = [&](int off) -> uint32_t {
            return (uint32_t)(data[off] | (data[off + 1] << 8) |
                (data[off + 2] << 16) | (data[off + 3] << 24));
            };

        auto Big = [&](int off, int n) -> BigInteger {
            std::vector<unsigned char> b(n + 1, 0);
            std::memcpy(b.data(), &data[off], n);
            // Convert little-endian / raw buffer to BigInteger string or bytes via GMP wrapper
            mpz_t temp_z;
            mpz_init(temp_z);
            mpz_import(temp_z, n + 1, -1, 1, 0, 0, b.data());
            BigInteger res;
            mpz_set(res.GetMpz(), temp_z);
            mpz_clear(temp_z);
            return res;
            };

        if (U32(0) != 0x44556677u) throw std::runtime_error("bad magic1");
        uint32_t fieldDataSize = U32(8);
        if (U32(12) != 0x00112233u) throw std::runtime_error("bad magic2");

        int pos = 20;
        unsigned char mustBe0 = data[pos];
        int sizeModulus = data[pos + 1];
        int sizeOrder = data[pos + 2];
        if (mustBe0 != 0) throw std::runtime_error("bad field data");

        int extDeg1 = (int)U32(pos + 3 + 0);
        int extDeg2 = (int)U32(pos + 3 + 4);
        pos = 20 + 3 + 36 + 1; // Skip 1 extra byte

        PubKey k;
        k.SizeModulus = sizeModulus;
        k.SizeOrder = sizeOrder;
        k.ExtDeg1 = extDeg1;
        k.ExtDeg2 = extDeg2;

        k.H1Bases.resize(sizeModulus);
        std::memcpy(k.H1Bases.data(), &data[pos], sizeModulus);
        pos += sizeModulus;

        k.Modulus = Big(pos, sizeModulus);
        pos += sizeModulus;
        Gf::SetP(k.Modulus); // Set global base field prime

        k.Order = Big(pos, sizeOrder);
        pos += sizeOrder;

        k.K3Minpoly.resize(extDeg1 + 1);
        std::memcpy(k.K3Minpoly.data(), &data[pos], extDeg1 + 1);
        pos += extDeg1 + 1;

        k.K6Minpoly.resize(extDeg2 + 1);
        std::memcpy(k.K6Minpoly.data(), &data[pos], extDeg2 + 1);
        pos += extDeg2 + 1;

        pos += sizeModulus * 2; // Skip padding/reserved section

        k.CurveA = Big(pos, sizeModulus);
        pos += sizeModulus;
        k.CurveB = Big(pos, sizeModulus);
        pos += sizeModulus;

        // Points and pairing values offset
        pos = (int)fieldDataSize + 12;
        for (int i = 0; i < sizeModulus; i++) {
            Fp3 x(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus), Big(pos + 2 * sizeModulus, sizeModulus));
            pos += sizeModulus * extDeg1;
            Fp3 y(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus), Big(pos + 2 * sizeModulus, sizeModulus));
            pos += sizeModulus * extDeg1;
            k.Points.push_back({ x, y });
        }

        Fp3 pv0(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus), Big(pos + 2 * sizeModulus, sizeModulus));
        pos += sizeModulus * extDeg1;
        Fp3 pv1(Big(pos, sizeModulus), Big(pos + sizeModulus, sizeModulus), Big(pos + 2 * sizeModulus, sizeModulus));
        k.PairingVal = { pv0, pv1 };

        return k;
    }
};

// -------------------------------------------------------------------------
// Product Key Calculation and Verification Engine
// -------------------------------------------------------------------------
class PKeyCalc {
private:
    static std::string Base64Encode(const unsigned char* bytes, size_t length) {
        static constexpr char kBase64Chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        std::string ret;
        int val = 0;
        int valb = -6;
        for (size_t i = 0; i < length; i++) {
            unsigned char c = bytes[i];
            val = (val << 8) + c;
            valb += 8;
            while (valb >= 0) {
                ret.push_back(kBase64Chars[(val >> valb) & 0x3F]);
                valb -= 6;
            }
        }
        if (valb > -6) ret.push_back(kBase64Chars[((val << 8) >> (valb + 8)) & 0x3F]);
        while (ret.size() % 4) ret.push_back('=');
        return ret;
    }

    static BigInteger ExtractM(const PubKey& k, const std::vector<unsigned char>& h1Coeffs) {
        int n = k.SizeModulus - 1;
        BigInteger acc = 0;
        for (int i = n - 1; i >= 0; i--) {
            int radix = k.H1Bases[1 + i] + 1;
            int digit = h1Coeffs[2 + i];
            if (digit >= radix) {
                throw std::runtime_error("Digit exceeds radix.");
            }
            acc = acc * radix + digit;
        }
        return acc;
    }

    static std::vector<unsigned char> ExtractMBytes(const PubKey& k, const std::vector<unsigned char>& h1Coeffs, int nBytes) {
        BigInteger mVal = ExtractM(k, h1Coeffs);
        mpz_t z;
        mpz_init(z);
        mpz_set(z, mVal.GetMpz());
        size_t count = 0;
        std::vector<unsigned char> raw(nBytes + 8, 0);
        mpz_export(raw.data(), &count, -1, 1, 0, 0, z);
        mpz_clear(z);

        std::vector<unsigned char> r(nBytes, 0);
        std::memcpy(r.data(), raw.data(), (std::min)(count, (size_t)nBytes));
        return r;
    }
    static bool TryCandidate(const PubKey& k, const std::vector<unsigned char>& enc, std::vector<unsigned char>& h1, std::vector<unsigned char>& uid) {
        std::vector<unsigned char> eb(enc.size() + 1, 0);
        std::memcpy(eb.data(), enc.data(), enc.size());

        mpz_t z;
        mpz_init(z);
        mpz_import(z, eb.size(), -1, 1, 0, 0, eb.data());
        BigInteger B;
        mpz_set(B.GetMpz(), z);
        mpz_clear(z);

        BigInteger xT = Gf::Mod(B);
        int keyByte = (int)((B / Gf::GetP()) & BigInteger(0xFF));

        std::optional<std::vector<BigInteger>> T_opt = Curve::LiftX(xT, k.CurveA, k.CurveB);
        if (!T_opt.has_value()) return false;
        std::vector<BigInteger> T = T_opt.value();

        h1.assign(15, 0);
        int b1 = k.H1Bases[1] + 1;
        int b2 = k.H1Bases[2] + 1;
        h1[0] = (unsigned char)keyByte;
        h1[1] = 1;
        h1[2] = (unsigned char)(keyByte % b1);
        h1[3] = (unsigned char)((keyByte / b1) % b2);

        Fp3 inv2 = Fp3(-2, 0, 0).Inverse();
        Fp3 inv2sq = inv2 * inv2;
        int np = (int)k.Points.size();

        // Shared Miller-ladder multi-pairing: the T-side doubling/addition
        // ladder is identical for every Q_i, so it's run once instead of
        // np times. Verified bit-identical to calling TatePairing::Pair
        // separately per point (see MultiPair in TatePairing.h).
        std::vector<Pt6> Qs(np);
        for (int i = 0; i < np; ++i) {
            Qs[i] = Pt6::From(Fp6::FromFp3(k.Points[i][0] * inv2),
                Fp6(Fp3::Zero(), k.Points[i][1] * inv2sq));
        }
        std::vector<Fp6> g = TatePairing::MultiPair(k.Order, xT, T[0], Qs, k.CurveA);

        std::vector<int> lo(np, 0), hi(np, 0);
        lo[0] = h1[1]; lo[1] = h1[2]; lo[2] = h1[3];
        for (int i = 0; i < 3; i++) hi[i] = lo[i];
        for (int i = 3; i < np; i++) hi[i] = k.H1Bases[i];

        Fp6 pv(k.PairingVal[0], k.PairingVal[1]);
        std::string note;
        // threads=0 -> draw from H1Search's shared global thread budget. When
        // this is the only entry doing real work (single-key validation) the
        // MITM parallelizes across all cores; when many entries run at once,
        // the budget keeps total threads near hardware_concurrency instead of
        // oversubscribing. Mirrors the C# reference's nested Parallel.For.
        std::vector<int> sol = H1Search::Solve(g, Fp6::One() / pv, lo, hi, note, 0);
        if (sol.empty()) return false;

        for (int i = 3; i < np; i++) {
            h1[4 + i - 3] = (unsigned char)sol[i];
        }

        uid = ExtractMBytes(k, h1, 8);
        return true;
    }

public:
    inline static double TLastPairMs = 0.0;
    inline static double TLastMitmMs = 0.0;

    // Inside class PKeyCalc (under public section):
    static bool TryParsedPubKey(const PubKey& k,
        const std::vector<unsigned char>& bEncryptArray,
        std::string& actPkeyConfig,
        std::vector<unsigned char>& h1Coeffs,
        std::vector<unsigned char>& uid) {
        actPkeyConfig.clear();
        h1Coeffs.clear();
        uid.clear();

        if (bEncryptArray.size() != 16) return false;

        std::vector<unsigned char> h1, u;
        if (!TryCandidate(k, bEncryptArray, h1, u)) return false;

        h1Coeffs = h1;
        uid = u;

        unsigned char cfg[12] = { 0 };
        std::memcpy(cfg, u.data(), (std::min)(u.size(), (size_t)8));
        actPkeyConfig = Base64Encode(cfg, 12);

        return true;
    }

    static bool TryPubKey(const std::vector<unsigned char>& publicKey,
        const std::vector<unsigned char>& bEncryptArray,
        std::string& actPkeyConfig,
        std::vector<unsigned char>& h1Coeffs,
        std::vector<unsigned char>& uid) {
        actPkeyConfig.clear();
        h1Coeffs.clear();
        uid.clear();

        if (publicKey.empty() || bEncryptArray.size() != 16) return false;

        PubKey k;
        try {
            k = PubKeyParser::Parse(publicKey);
        }
        catch (const BaseFieldMismatchException&) {
            throw;
        }
        catch (...) {
            return false;
        }

        std::vector<unsigned char> h1, u;
        if (!TryCandidate(k, bEncryptArray, h1, u)) return false;

        h1Coeffs = h1;
        uid = u;

        unsigned char cfg[12] = { 0 };
        std::memcpy(cfg, u.data(), (std::min)(u.size(), (size_t)8));
        actPkeyConfig = Base64Encode(cfg, 12);

        int ev = u[0] | ((u[1] | ((u[2] | ((u[3] & 0x7F) << 8)) << 8)) << 8);
        return true;
    }
};
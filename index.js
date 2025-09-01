// solve_robust_newton.js
// Usage: node solve_robust_newton.js path/to/input.json
// Robust solver for c = P(0) using Newton's interpolation (no Lagrange).
// - Tries all C(n, k) subsets
// - Picks the polynomial (degree k-1) that fits all n points exactly
// - Uses BigInt + rational arithmetic for exactness

const fs = require("fs");

// ---------- utils ----------
function die(msg) { console.error(msg); process.exit(1); }

// gcd for BigInt
function bigIntGcd(a, b) {
    a = a < 0n ? -a : a;
    b = b < 0n ? -b : b;
    while (b !== 0n) { const t = a % b; a = b; b = t; }
    return a;
}

function normSign(num, den) {
    if (den < 0n) { num = -num; den = -den; }
    return [num, den];
}

// Rational as reduced fraction of BigInts
function makeRat(num, den = 1n) {
    num = BigInt(num); den = BigInt(den);
    [num, den] = normSign(num, den);
    const g = bigIntGcd(num, den);
    return { num: num / g, den: den / g };
}
function rAdd(a, b) { return makeRat(a.num * b.den + b.num * a.den, a.den * b.den); }
function rSub(a, b) { return makeRat(a.num * b.den - b.num * a.den, a.den * b.den); }
function rMul(a, b) { return makeRat(a.num * b.num, a.den * b.den); }
function rDiv(a, b) {
    if (b.num === 0n) throw new Error("Division by zero");
    return makeRat(a.num * b.den, a.den * b.num);
}
function rEqInt(a, nBig) { return a.den === 1n && a.num === nBig; }

// base-N string -> BigInt (supports up to base 36)
function toBigIntFromBase(valueStr, baseStr) {
    const base = BigInt(parseInt(baseStr, 10));
    if (base < 2n || base > 36n) throw new Error(`Unsupported base ${baseStr}`);
    const digits = "0123456789abcdefghijklmnopqrstuvwxyz";
    const map = new Map(digits.split("").map((ch, i) => [ch, BigInt(i)]));
    let acc = 0n;
    for (const raw of valueStr.toLowerCase()) {
        if (!map.has(raw)) throw new Error(`Invalid digit '${raw}'`);
        const d = map.get(raw);
        if (d >= base) throw new Error(`Digit '${raw}' not valid for base ${baseStr}`);
        acc = acc * base + d;
    }
    return acc;
}

// Build Newton coefficients (divided differences) for given points
// returns array of rationals: coeffs[0..k-1]
function newtonCoeffs(xs, ysInt) {
    const k = xs.length;
    // as rationals
    const col0 = ysInt.map(y => makeRat(y, 1n));
    const coeffs = col0.slice();
    for (let d = 1; d < k; d++) {
        for (let i = 0; i < k - d; i++) {
            const num = rSub(coeffs[i + 1], coeffs[i]);
            const den = makeRat(xs[i + d] - xs[i], 1n);
            coeffs[i] = rDiv(num, den);
        }
    }
    // After in-place building, coeffs[0..k-1] are the Newton coefficients
    return coeffs.slice(0, k);
}

// Evaluate Newton form at X (as rational)
function newtonEvalAt(coeffs, xs, X) {
    let result = coeffs[0];
    let prod = makeRat(1n, 1n);
    for (let i = 1; i < coeffs.length; i++) {
        prod = rMul(prod, makeRat(X - xs[i - 1], 1n));
        result = rAdd(result, rMul(coeffs[i], prod));
    }
    return result;
}

// Generate all k-size subsets of indices [0..n-1]
function* kSubsets(n, k) {
    const idx = Array.from({ length: k }, (_, i) => i);
    while (true) {
        yield idx.slice();
        let i = k - 1;
        while (i >= 0 && idx[i] === i + n - k) i--;
        if (i < 0) break;
        idx[i]++;
        for (let j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
    }
}

// ---------- main ----------
const file = process.argv[2] || "sample1.json";
let obj;
try { obj = JSON.parse(fs.readFileSync(file, "utf8")); }
catch (e) { die("Failed to read/parse JSON: " + e.message); }

const k = obj?.keys?.k;
const n = obj?.keys?.n;
if (!Number.isInteger(k) || !Number.isInteger(n)) die("Invalid keys.n/keys.k");

const allPoints = Object.keys(obj)
    .filter(k => k !== "keys")
    .map(kStr => {
        const x = BigInt(parseInt(kStr, 10));
        const { base, value } = obj[kStr];
        const y = toBigIntFromBase(value, base);
        return [x, y];
    })
    .sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0));

if (allPoints.length !== n) {
    // Not fatal, but warn
    console.error(`Warning: keys.n=${n} but found ${allPoints.length} data points`);
}

const xsAll = allPoints.map(p => p[0]);
const ysAll = allPoints.map(p => p[1]);

let cCandidate = null; // BigInt
let foundCount = 0;

for (const comb of kSubsets(allPoints.length, k)) {
    const xs = comb.map(i => xsAll[i]);
    const ys = comb.map(i => ysAll[i]);

    let coeffs;
    try {
        coeffs = newtonCoeffs(xs, ys);
    } catch (e) {
        // singular configuration; skip
        continue;
    }
    // Evaluate at all points; must equal integer y exactly
    let ok = true;
    for (let i = 0; i < xsAll.length; i++) {
        const val = newtonEvalAt(coeffs, xs, xsAll[i]);
        if (!rEqInt(val, ysAll[i])) { ok = false; break; }
    }
    if (!ok) continue;

    // Good polynomial; take c = P(0)
    const cRat = newtonEvalAt(coeffs, xs, 0n);
    if (cRat.den !== 1n) {
        // For exact integer polynomials, this shouldn't happen
        continue;
    }
    const cVal = cRat.num;

    if (cCandidate === null) {
        cCandidate = cVal;
    } else if (cCandidate !== cVal) {
        die(`Ambiguity detected: multiple valid polynomials with different c (e.g., ${cCandidate} vs ${cVal})`);
    }
    foundCount++;
}

// If no subset fits all points, fall back to using the first k points (still returns something)
if (cCandidate === null) {
    const xs = xsAll.slice(0, k);
    const ys = ysAll.slice(0, k);
    const coeffs = newtonCoeffs(xs, ys);
    const cRat = newtonEvalAt(coeffs, xs, 0n);
    if (cRat.den !== 1n) die(`Non-integer c from fallback: ${cRat.num}/${cRat.den}`);
    console.error("Note: no subset fit all points exactly; using first k points.");
    console.log(String(cRat.num));
} else {
    // Optionally report how many subsets matched
    // console.error(`Matched subsets: ${foundCount}`);
    console.log(String(cCandidate));
}

// Downloads the processed gait files of Fukuchi, Fukuchi & Duarte (2017) into data/fukuchi2017/.
// Source: Figshare article 4543435 (doi 10.6084/m9.figshare.4543435), CC BY 4.0. Run: node tools/fetch-fukuchi.js
// With --raw it also downloads, for subjects 1-28, the raw marker and force files of the three
// running speeds and the static trial (RBDSnnnrunT25markers.txt, ...forces.txt, RBDSnnnstatic.txt; about 330 MB).
'use strict';
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const ARTICLE = 'https://api.figshare.com/v2/articles/4543435';
const DEST = path.join(__dirname, '..', 'data', 'fukuchi2017');
const RAW = process.argv.includes('--raw');
const CONCURRENCY = 4, RETRIES = 4;

const md5 = buf => crypto.createHash('md5').update(buf).digest('hex');
const sleep = ms => new Promise(r => setTimeout(r, ms));

function rawNames() {
  const names = new Set();
  for (let s = 1; s <= 28; s++) {
    const id = 'RBDS' + String(s).padStart(3, '0');
    names.add(id + 'static.txt');
    for (const sp of ['25', '35', '45']) { names.add(`${id}runT${sp}markers.txt`); names.add(`${id}runT${sp}forces.txt`); }
  }
  return names;
}

// returns 'present' | 'downloaded'; throws after all retries fail
async function fetchOne(f) {
  const dest = path.join(DEST, f.name);
  if (fs.existsSync(dest)) {
    const buf = fs.readFileSync(dest);
    if (buf.length === f.size && md5(buf) === f.computed_md5) return { status: 'present', bytes: buf.length };
  }
  let lastErr;
  for (let attempt = 1; attempt <= RETRIES; attempt++) {
    try {
      const r = await fetch(`https://ndownloader.figshare.com/files/${f.id}`, { redirect: 'follow' });
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const buf = Buffer.from(await r.arrayBuffer());
      if (buf.length !== f.size) throw new Error(`size ${buf.length}, expected ${f.size}`);
      if (md5(buf) !== f.computed_md5) throw new Error('md5 mismatch');
      fs.writeFileSync(dest, buf);
      return { status: 'downloaded', bytes: buf.length };
    } catch (e) {
      lastErr = e;
      if (attempt < RETRIES) await sleep(1000 * attempt);
    }
  }
  throw lastErr;
}

async function main() {
  fs.mkdirSync(DEST, { recursive: true });
  const res = await fetch(ARTICLE);
  if (!res.ok) throw new Error(`listing: HTTP ${res.status}`);
  const files = (await res.json()).files;
  const want = files.filter(f => f.name.endsWith('processed.txt'));
  const infos = files.filter(f => f.name === 'RBDSinfo.txt').sort((a, b) => b.id - a.id);
  if (infos.length) want.push(infos[0]);
  else throw new Error('RBDSinfo.txt not in the listing');
  if (RAW) {
    const names = rawNames(), found = new Set();
    for (const f of files) if (names.has(f.name) && !found.has(f.name)) { want.push(f); found.add(f.name); }
    const missing = [...names].filter(n => !found.has(n));
    if (missing.length) console.warn(`not in the listing (${missing.length}): ${missing.join(', ')}`);
  }

  let present = 0, downloaded = 0, failed = 0, bytes = 0, next = 0;
  async function worker() {
    while (next < want.length) {
      const f = want[next++];
      try {
        const r = await fetchOne(f);
        bytes += r.bytes;
        if (r.status === 'present') present++;
        else { downloaded++; console.log(`downloaded ${f.name} (${r.bytes} bytes)`); }
      } catch (e) {
        failed++;
        console.error(`FAILED ${f.name}: ${e.message}`);
      }
    }
  }
  await Promise.all(Array.from({ length: RAW ? CONCURRENCY : 1 }, worker));
  console.log(`${want.length} files: ${present} already present and verified, ${downloaded} downloaded, ${failed} failed; ${(bytes / 1e6).toFixed(1)} MB verified in total`);
  if (failed) process.exit(1);
}

main().catch(e => { console.error(e.message); process.exit(1); });

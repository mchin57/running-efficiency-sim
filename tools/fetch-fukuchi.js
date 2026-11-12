// Downloads the processed gait files of Fukuchi, Fukuchi & Duarte (2017) into data/fukuchi2017/.
// Source: Figshare article 4543435 (doi 10.6084/m9.figshare.4543435), CC BY 4.0. Run: node tools/fetch-fukuchi.js
'use strict';
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const ARTICLE = 'https://api.figshare.com/v2/articles/4543435';
const DEST = path.join(__dirname, '..', 'data', 'fukuchi2017');

const md5 = buf => crypto.createHash('md5').update(buf).digest('hex');

async function main() {
  fs.mkdirSync(DEST, { recursive: true });
  const res = await fetch(ARTICLE);
  if (!res.ok) throw new Error(`listing: HTTP ${res.status}`);
  const files = (await res.json()).files;
  const want = files.filter(f => f.name.endsWith('processed.txt'));
  const infos = files.filter(f => f.name === 'RBDSinfo.txt').sort((a, b) => b.id - a.id);
  if (infos.length) want.push(infos[0]);
  else throw new Error('RBDSinfo.txt not in the listing');

  let present = 0, downloaded = 0, failed = 0;
  for (const f of want) {
    const dest = path.join(DEST, f.name);
    try {
      if (fs.existsSync(dest)) {
        const buf = fs.readFileSync(dest);
        if (buf.length === f.size && md5(buf) === f.computed_md5) { present++; continue; }
      }
      const r = await fetch(`https://ndownloader.figshare.com/files/${f.id}`, { redirect: 'follow' });
      if (!r.ok) throw new Error(`HTTP ${r.status}`);
      const buf = Buffer.from(await r.arrayBuffer());
      if (buf.length !== f.size) throw new Error(`size ${buf.length}, expected ${f.size}`);
      if (md5(buf) !== f.computed_md5) throw new Error('md5 mismatch');
      fs.writeFileSync(dest, buf);
      downloaded++;
      console.log(`downloaded ${f.name} (${buf.length} bytes)`);
    } catch (e) {
      failed++;
      console.error(`FAILED ${f.name}: ${e.message}`);
    }
  }
  console.log(`${want.length} files: ${present} already present and verified, ${downloaded} downloaded, ${failed} failed`);
  if (failed) process.exit(1);
}

main().catch(e => { console.error(e.message); process.exit(1); });

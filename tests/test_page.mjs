/* Offline tests for the website.

   The whole <script> from site/index.html runs against a stubbed DOM and a
   stubbed fetch serving synthetic genes, so there are no fixture files and no
   network. Run: node tests/test_page.mjs                                   */

import fs from 'fs';
import path from 'path';

/* ---- synthetic data ----------------------------------------------------
   PLUSG on chr7 (+): coding exons 1-3 at 1000-1099 / 2000-2099 / 3000-3099.
   Exon 1 has no primers, exon 2 is OK, exon 3 is CHECK (a SNP under F).
   MINUSG on chr8 (-) with the same coordinates and reversed ranks.
   REG1 is a custom region on chr7:5000-5300.                            */
function pair(chrom, start, end, extra={}){
  const amp = 'ACGTACGTACGTACGTACGT' + 'acgt'.repeat(10) + 'GGGTTTCCCAAAGGGTTTCC';
  return {rank:1, f_seq:'ACGTACGTACGTACGTACGT', r_seq:'GGAAACCCTTTGGGAAACCC',
    f_tm:59.1, r_tm:58.7, f_tm_ucsc:63.2, r_tm_ucsc:62.8, d_tm:0.4, f_gc:50, r_gc:50,
    f_end_gc:2, r_end_gc:2, f_len:20, r_len:20, product:amp.length, product_tm:81.4,
    penalty:0.62, placement:'exonic', target_cov:80, chrom, start, end, amplicon:amp,
    f_snps:[], r_snps:[], amp_snps:[], f_repeat:0, r_repeat:0, repeat_names:[],
    ispcr:{checked:true,on_target:true,off_target:[]}, ...extra};
}
function ampRec(gene, n, chrom, start, status, extra={}){
  const warnings = status==='OK' ? [] : ['F primer: SNP rs123 (MAF 0.12) 8 nt from 3\' end <script>'];
  return {label:String(n), name:`${gene}_qE${n}`, exons:[n], target_start:start, target_end:start+99,
    target_len:100, exon_len:100, min_cov:40, relaxed:false, tier:status==='OK'?'clean':'snp-5p',
    placement:'exonic', segdups:[], status, warnings,
    notes: status==='OK' ? [] : ['1 common SNP(s) between the primers (not under them): rs456'],
    pairs:[pair(chrom, start+10, start+89, status==='OK'?{}:{f_snps:[{rs:'rs123',maf:0.12,pos3:8}]}),
           {...pair(chrom, start+5, start+84), rank:2, penalty:0.9}], ...extra};
}
function mkGene(name, chrom, strand){
  const base = [{start:1000,end:1099},{start:2000,end:2099},{start:3000,end:3099}];
  const order = strand===1 ? base : [...base].reverse();
  const exons = order.map((e,i)=>({n:i+1, start:e.start, end:e.end, cds_start:e.start, cds_end:e.end,
                                     coding:true, len:100}));
  exons.push({n:4, start:strand===1?4000:500, end:strand===1?4099:599, cds_start:null, cds_end:null, coding:false, len:100});
  const at = n => exons.find(e=>e.n===n).start;
  return {gene:name, kind:'gene', gene_id:'ENSG0', transcript:`ENST_${name}.1`, refseq:`NM_${name}.1`,
    selection:'MANE Select', assembly:'GRCh38', chrom, strand, exon_count:4,
    annotated:{snps:true,repeats:true,segdups:true}, exons, n_snps:3,
    amplicons:[ampRec(name, 2, chrom, at(2), 'OK'), ampRec(name, 3, chrom, at(3), 'CHECK')],
    failures:[{label:'1', name:`${name}_qE1`, exons:[1], target_start:at(1), target_end:at(1)+99,
      target_len:100, exon_len:100, min_cov:40, relaxed:false, tier:'', placement:'', segdups:[],
      reason:'left: considered 10, GC content failed 10'}]};
}
const PLUSG = mkGene('PLUSG', '7', 1);
const MINUSG = mkGene('MINUSG', '8', -1);
const REG1 = {gene:'REG1', kind:'region', gene_id:'', transcript:'', refseq:'', selection:'custom region',
  assembly:'GRCh38', chrom:'7', strand:1, exon_count:1, annotated:{snps:true,repeats:true,segdups:true},
  exons:[{n:1, start:5000, end:5300, cds_start:5000, cds_end:5300, coding:true, len:301}], n_snps:0,
  amplicons:[ampRec('REG1', 1, '7', 5000, 'OK', {label:'', name:'REG1_q'})], failures:[]};

const entry = g => ({gene:g.gene, kind:g.kind, transcript:g.transcript, refseq:g.refseq, chrom:g.chrom,
  start:Math.min(...g.exons.map(e=>e.start)), end:Math.max(...g.exons.map(e=>e.end)),
  amplicons:g.amplicons.length, ok:1, check:0, failures:g.failures.length});
const INDEX = {
  assembly:'GRCh38', ucsc_db:'hg38', generated:'2026-11-12', version:'1.0.0',
  params:{target:'coding', min_overlap:40, flank:100, product_min:50, product_max:150, product_opt:80,
    tm_min:58, tm_opt:59, tm_max:60, size_min:18, size_opt:20, size_max:28, gc_min:30, gc_max:80,
    max_end_gc:2, max_poly_x:3, max_tm_diff:2, snp_maf:0.01, snp_3prime:5, salt_mono:50, salt_div:1.5,
    dntp:0.6, dna_conc:50, num_return:3},
  ispcr:true, tiers:{}, errors:[{gene:'BROKEN', error:'lookup failed'}],
  genes:[entry(PLUSG), entry(MINUSG), entry(REG1)],
};
const FILES = {'data/index.json':INDEX, 'data/genes/PLUSG.json':PLUSG,
               'data/genes/MINUSG.json':MINUSG, 'data/genes/REG1.json':REG1};

/* ---- stub DOM + fetch, then run the real page script ---- */
const els = {};
const el = id => (els[id] = els[id] || {id, innerHTML:'', textContent:'', value:'',
  classList:{add(){},remove(){}}, addEventListener(){}, dataset:{}, onclick:null});
globalThis.document = {querySelector: s => el(s.replace('#','')), querySelectorAll: () => []};
globalThis.location = {hash:'', pathname:'/qpcrprimer-batch/', search:''};
globalThis.history = {pushState(){}};
globalThis.window = {addEventListener(){}};
Object.defineProperty(globalThis, 'navigator',
  {value:{clipboard:{writeText:()=>Promise.resolve()}}, configurable:true});
let fetches = [];
globalThis.fetch = async u => {
  fetches.push(u);
  if(u in FILES) return {json: async () => JSON.parse(JSON.stringify(FILES[u]))};
  throw new Error('404 ' + u);
};

const page = fs.readFileSync(path.join(import.meta.dirname, '..', 'site', 'index.html'), 'utf8');
const script = page.match(/<script>([\s\S]*)<\/script>/);
if(!script){ console.error('no <script> block in site/index.html'); process.exit(1); }
const api = new Function(script[1] + `
  return {init, route, showRegion, parseCoords, exonsInRange, submitQuery, renderGene,
          geneTSV, geneOrder, orderPair, loadGene,
          view:()=>document.querySelector('#view').innerHTML};`)();

let pass = 0, fail = 0;
const ok = (c, m) => { if(c) pass++; else { fail++; console.log('  FAIL: ' + m); } };
const view = () => api.view();
const has  = t => view().includes(t);
const go = async h => { location.hash = '#' + h; await api.route(); };

await api.init();

/* ---- init ---- */
ok(els['meta'].innerHTML.includes('2 genes') && els['meta'].innerHTML.includes('1 region'),
   'header counts genes and custom regions separately');
ok(els['genelist'].innerHTML.includes('PLUSG') && els['genelist'].innerHTML.includes('custom region'),
   'datalist lists genes and regions');
ok(els['note'].innerHTML.includes('Primer Express'), 'method note names the rules');
ok(els['note'].innerHTML.includes('58–60'), 'and the Tm window');
ok(els['note'].innerHTML.includes('BROKEN'), 'genes that failed are listed');
ok(!els['note'].innerHTML.includes('undefined'), 'no undefined in the method note');

/* ---- gene view ---- */
await go('PLUSG');
ok(has('class="ex coding fail sel" data-n="1"'), 'exon 1 (no primers) is selected by default and marked failed');
ok(has('class="ex coding ok" data-n="2"'), 'exon 2 is OK');
ok(has('class="ex coding check" data-n="3"'), 'exon 3 is Check');
ok(has('title="Exon 4: non-coding"'), 'non-coding exon is shown but not clickable');
ok(has('PLUSG_qE1') && has('No primers') && has('GC content failed'), 'a failed exon says why');

await go('PLUSG/2');
ok(has('PLUSG_qE2F') && has('PLUSG_qE2R'), 'primer names are GENE_qE<n>F/R');
ok(has('ACGTACGTACGTACGTACGT') && has('GGAAACCCTTTGGGAAACCC'), 'both primers are shown');
ok(has('81.4 °C'), 'product Tm (melt peak) is shown');
ok(has('0.62'), 'penalty is shown');
ok(has('inside the exon'), 'placement is stated');
ok(has('2/5'), "3' G/C count is shown");
ok(!has('class="hot"'), 'nothing is highlighted when every value is inside the rules');
ok(has('1 alternative pair'), 'alternatives are offered');
ok(has('hgPcr?') && has('wp_f=ACGTACGTACGTACGTACGT'), 'UCSC in-silico PCR link carries the pair');
ok(!has('undefined') && !has('NaN'), 'no undefined or NaN in the exon view');

await go('PLUSG/3');
ok(has('pill check'), 'a Check amplicon carries the Check pill');
ok(has('&lt;script&gt;') && !has('<script>'), 'warnings are escaped');
ok(has('rs456'), 'notes are shown');
ok(has('rs123') && has('8 nt from 3'), 'the SNP under a primer is shown under the sequence');

FILES['data/genes/MINUSG.json'].amplicons[0].pairs[0].f_tm = 57.2;   // as a relaxed pair can be
await go('MINUSG/2');
ok(has('<td class="hot">57.2</td>'), 'a Tm outside the window is highlighted');
ok(has('MINUSG_qE2') && has('chr<b>8</b> (−)'), 'minus-strand gene renders');

await go('REG1');
ok(has('custom region') && has('chr7:5,000-5,300'), 'a custom region renders with its span');
ok(has('REG1_qF') && !has('class="ex coding'), 'with its primers and no exon chips');
ok(has('inside the region'), 'placement speaks of the region');

await go('NOTAGENE');
ok(has('not in this build'), 'an unknown symbol is reported');
await go('PLUSH');
ok(has('Did you mean') && has('>PLUSG<'), 'a one-letter typo suggests the real symbol');

/* ---- coordinates ---- */
const c1 = api.parseCoords('chr7:1,995-2,015');
ok(c1 && c1.chrom==='7' && c1.start===1995 && c1.end===2015 && !c1.point, 'range with commas parses');
const c2 = api.parseCoords('X:300-100');
ok(c2 && c2.chrom==='X' && c2.start===100 && c2.end===300, 'reversed range is ordered');
ok(api.parseCoords('chr7:2050').point, 'a single position is a point');
ok(api.parseCoords('chrM:5').chrom==='MT', 'chrM is MT');
ok(api.parseCoords('PTEN')===null && api.parseCoords('chr99:1')===null && api.parseCoords('chr7:a-b')===null,
   'non-coordinates are refused');

fetches = [];
await go('q=chr7:1995-2015');
ok(has('PLUSG_qE2') && !has('PLUSG_qE3'), 'a range inside exon 2 lists exon 2 only');
ok(!fetches.some(u=>u.includes('MINUSG')), 'genes on another chromosome are not fetched');

await go('q=chr7:900-3200');
ok(has('PLUSG_qE1') && has('PLUSG_qE2') && has('PLUSG_qE3'), 'a range over the gene lists every coding exon');
ok(has('<b>3</b> assays, 1 OK'), 'and counts the assays');
ok(has('No primers'), 'a failed exon in the range is listed, not hidden');

await go('q=chr7:2500');
ok(has('intronic in PLUSG') && has('PLUSG_qE2') && has('PLUSG_qE3'), 'an intronic point shows the exons on either side');

await go('q=chr8:2050');
ok(has('MINUSG_qE2'), 'minus-strand gene found by position');

await go('q=chr7:5100-5200');
ok(has('REG1_q') && has('REG1'), 'a custom region is found by position');
ok(has('#REG1"'), 'and links to the region view');

await go('q=chr22:1-100');
ok(has('No gene or custom region'), 'a range outside the build says so');
ok(has('regions.txt'), 'and says how to add it');

/* ---- routing ---- */
location.hash = '';
api.submitQuery('plusg');
ok(location.hash === 'PLUSG', `gene symbol routes as a gene (got "${location.hash}")`);
api.submitQuery('chr7:1,000 - 2,000');
ok(location.hash === 'q=' + encodeURIComponent('chr7:1,000-2,000'), `region routes as a search (got "${location.hash}")`);
api.submitQuery('reg1');
ok(location.hash === 'REG1', 'a custom region name routes to it');

/* ---- exports ---- */
const g = await api.loadGene('PLUSG');
const rows = api.orderPair(g, g.amplicons[0], g.amplicons[0].pairs[0]).split('\n');
ok(rows.length === 2, 'one order line per oligo');
const cols = rows[0].split('\t');
ok(cols.length === 7, `seven columns (got ${cols.length})`);
ok(cols[0] === 'PLUSG-qE2-F' && rows[1].split('\t')[0] === 'PLUSG-qE2-R', `names are GENE-qE<n>-F/R (got "${cols[0]}")`);
ok(cols[1] === 'ACGTACGTACGTACGTACGT' && cols[2] === '20' && cols[3] === '', 'sequence, length, empty control cell');
ok(cols[4] === '40 nmol (10-99 bases)' && cols[5] === 'RP-HPLC' && cols[6] === 'H2O-100 µM', 'fixed order parameters');
ok(api.geneOrder(g).split('\n').length === 4, 'gene order rows: 2 amplicons x 2 oligos');
const r = await api.loadGene('REG1');
ok(api.orderPair(r, r.amplicons[0], r.amplicons[0].pairs[0]).startsWith('REG1-q-F\t'), 'a region is named REGION-q-F');
const tsv = api.geneTSV(g).split('\n');
ok(tsv[0].startsWith('Name\tSequence') && tsv.length === 5, 'TSV: header + 2 lines per amplicon');
ok(tsv[1].startsWith('PLUSG_qE2F\tACGT'), 'TSV lines carry name and sequence');

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);

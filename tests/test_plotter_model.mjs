import assert from 'node:assert/strict';
import {makeTraces,maskForLog,spectrumSteps,records,csvText,filteredSamples,axisRange} from '../plotter/static/model.mjs';

assert.deepEqual(axisRange(null,null),{range:[null,null],autorange:true});
assert.deepEqual(axisRange(10,null),{range:[10,null],autorange:'max'});
assert.deepEqual(axisRange(null,100),{range:[null,100],autorange:'min'});
assert.deepEqual(axisRange(1,100,true),{range:[0,2],autorange:false});
assert.deepEqual(axisRange(null,1000,true),{range:[null,3],autorange:'min'});
assert.throws(()=>axisRange(0,null,true));
assert.throws(()=>axisRange(10,10));
assert.throws(()=>axisRange(100,10));
assert.throws(()=>axisRange(NaN,null));

const samples=[{time:100,label:'100 seconds',zone:1,set:6},
               {time:200,label:'200 seconds',zone:1,set:6}];
const cell={cell:'6000',families:{activity:{unit:'Bq/cm³',samples,
  series:[{name:'H3',values:[0,20]},{name:'TOTAL',values:[10,30]}]},
  gamma:{unit:'photons/(cm³·s)',samples,bins:[{lower:0,upper:.01},{lower:.01,upper:.02}],
    values:[[5,10],[0,2]],totals:[15,2],series:[]}}};
const config={family:'activity',mode:'evolution',timeUnit:'days',zone:'all',timeSet:'all',
  minimum:null,maximum:null,series:new Set(['H3']),times:new Set(['100']),bins:new Set(['TOTAL'])};
const trace=makeTraces([cell],config)[0];
assert.ok(!trace.name.includes(' · Z'));
assert.deepEqual(trace.points.map(p=>p.y),[0,20]);
assert.equal(trace.points[0].x,100/86400);
const masked=maskForLog([trace],true,true);
assert.equal(masked.masked,1);assert.equal(masked.visible,1);
assert.equal(masked.series[0].points[0].y,null);
assert.equal(trace.points[0].y,0); // Log scale does not alter export/source values.
assert.equal(records([trace],'activity')[0].value,0);
assert.equal(filteredSamples(cell.families.activity,{...config,minimum:150}).length,1);
const spectral=makeTraces([cell],{...config,family:'gamma',mode:'spectrum'})[0];
assert.deepEqual(spectral.points.map(p=>p.x),[.01,.02]);
assert.deepEqual(spectral.points.map(p=>p.y),[5,10]);
assert.equal(records([spectral],'gamma')[0].lower_energy_MeV,0);
assert.equal(records([spectral],'gamma')[1].upper_energy_MeV,.02);
const steps=spectrumSteps([{x:.01,lower:0,y:null},{x:.02,lower:.01,y:2}],false);
assert.deepEqual(steps.map(p=>[p.x,p.y]),[[0,null],[.01,null],[.01,2],[.02,2]]);
assert.ok(spectrumSteps(spectral.points,true).every(p=>p.x>0));
const temporal=makeTraces([cell],{...config,family:'gamma',bins:new Set(['.01','0.02'])});
assert.equal(temporal.length,1); // Keys use their exact numeric string, not approximate matches.
assert.deepEqual(temporal[0].points.map(p=>p.y),[10,2]);
const duplicate=structuredClone(cell);
duplicate.families.activity.samples=[samples[0],{...samples[0],set:7}];
assert.equal(makeTraces([duplicate],config).length,2); // Distinct restarts cannot merge.
duplicate.families.activity.samples=[samples[0],{...samples[0],zone:2}];
assert.equal(makeTraces([duplicate],config).length,2);
assert.ok(makeTraces([duplicate],config).every(t=>t.name.includes(' · Z')));
const csv=csvText([{cell:'=unsafe',value:0,text:'one, "two"',negative:-1}]);
assert.ok(csv.includes("'=unsafe"));assert.ok(csv.includes('"one, ""two"""'));assert.ok(csv.endsWith(',-1'));
console.log('Plotter model: filters, units, bins, duplicates, log scales and CSV passed.');

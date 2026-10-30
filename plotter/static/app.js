import {units,labels,filteredSamples,rankedSeries,makeTraces,maskForLog,spectrumSteps,records,csvText,axisRange} from './model.mjs';

const $=id=>document.getElementById(id);
const state={cells:[],chosen:new Set(),cache:new Map(),family:'activity',mode:'spectrum',
  series:new Set(['TOTAL']),times:new Set(),bins:new Set(['TOTAL']),zone:'all',timeSet:'all',
  timeUnit:'seconds',minimum:null,maximum:null,view:'chart',raw:[],request:0,curveColors:new Map(),limits:new Map(),fontSize:14,legendSize:14,
  axes:{evolution:{x:true,y:true},spectrum:{x:false,y:true}}};
const colors=['#087e8b','#e18a39','#6974b1','#bd5b73','#329c77','#8d6a43','#496c9a','#9a69a6','#6d8c3e','#ce6844'];
const dashes=['solid','dash','dot','dashdot','longdash','longdashdot'];
const friendly=name=>name==='TOTAL'?'Printed total':name;
const number=value=>Number(value).toLocaleString('en-US',{maximumSignificantDigits:4});
const scientific=value=>Number(value).toExponential(3);
const escape=value=>String(value).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;');
const spectrum=()=>state.family==='gamma'&&state.mode==='spectrum';
const axis=()=>state.axes[spectrum()?'spectrum':'evolution'];
const limitKey=()=>`${state.family}:${spectrum()}`;
const limits=()=>state.limits.get(limitKey())??{xmin:null,xmax:null,ymin:null,ymax:null};
const xFactor=()=>spectrum()?1:units[state.timeUnit].factor;
function drawAxisLimits() {
  const values=limits();
  for(const key of ['xmin','xmax','ymin','ymax'])$(key).value=values[key]===null?'':values[key]/(key.startsWith('x')?xFactor():1);
  $('axis-x-label').textContent=`X [${spectrum()?'MeV':units[state.timeUnit].label}]`;
  $('axis-limit-error').hidden=true;
}
const loaded=()=>state.cells.filter(c=>state.chosen.has(c.id)&&state.cache.has(c.id)).map(c=>state.cache.get(c.id));
const config=()=>({family:state.family,mode:state.mode,series:state.series,times:state.times,bins:state.bins,
  zone:state.zone,timeSet:state.timeSet,timeUnit:state.timeUnit,minimum:state.minimum,maximum:state.maximum});

let viewerToken=null,viewerConnection=null;
function connectViewer(token=viewerToken) {
  if(!token)return;
  if(viewerConnection&&viewerToken===token)return;
  viewerConnection?.close();viewerToken=token;
  viewerConnection=new EventSource(`/api/viewer/${token}`);
}
window.addEventListener('pagehide',()=>{viewerConnection?.close();viewerConnection=null;});
window.addEventListener('pageshow',event=>{
  if(event.persisted)connectViewer();
});

async function api(path,options={}) {
  const response=await fetch(path,{cache:'no-store',...options});
  const data=await response.json();
  if (!response.ok) throw new Error(data.error||`HTTP error ${response.status}`);
  return data;
}
const post=(path,data)=>api(path,{method:'POST',headers:{'Content-Type':'application/json','X-APYPA-Token':viewerToken},body:JSON.stringify(data)});
function busyProject(busy) {
  for(const id of ['pick-directory','load-directory','export-excel','project-directory'])$(id).disabled=busy;
}
async function loadProject(directory) {
  if(!directory.trim())throw new Error('Enter a directory.');
  await post('/api/project',{directory});
  state.cache.clear();state.chosen.clear();state.series=new Set(['TOTAL']);
  $('excel-download').hidden=true;await discover();
  $('project-status').textContent='Directory loaded.';
}
$('load-directory').onclick=async()=>{
  busyProject(true);error([]);
  try {await loadProject($('project-directory').value);}catch(e){error([e.message]);}
  finally {busyProject(false);}
};
$('project-directory').addEventListener('keydown',event=>{if(event.key==='Enter')$('load-directory').click();});
$('pick-directory').onclick=async()=>{
  busyProject(true);error([]);$('project-status').textContent='Choose a directory in the system dialog…';
  try {const result=await post('/api/project/pick',{});if(result.directory)await loadProject(result.directory);else $('project-status').textContent='Selection canceled.';}
  catch(e){error([e.message]);$('project-status').textContent='You can enter the path and click Load.';}
  finally {busyProject(false);}
};
$('export-excel').onclick=async()=>{
  busyProject(true);error([]);$('excel-download').hidden=true;
  $('project-status').textContent='Generating Excel workbooks for the selected cells…';
  try {
    const result=await post('/api/export',{cells:[...state.chosen],overwrite:$('overwrite-excel').checked});
    $('project-status').textContent=`${result.files.length} Excel ${result.files.length===1?'workbook':'workbooks'} saved in ${result.directory}`;
    if(result.errors.length)error(result.errors);
    if(result.download){$('excel-download').href=result.download;$('excel-download').hidden=false;}
  } catch(e){error([e.message]);$('project-status').textContent='The export did not complete.';}
  finally {busyProject(false);}
};
function error(messages) {
  $('errors').hidden=!messages.length;
  $('errors').textContent=messages.join(' · ');
}
function option(container, value, text, checked, callback, secondary='') {
  const label=document.createElement('label');label.className=`option${checked?' selected':''}`;
  const input=document.createElement('input');input.type='checkbox';input.value=value;input.checked=checked;
  input.addEventListener('change',()=>{label.classList.toggle('selected',input.checked);callback(input.checked);});
  const title=document.createElement('span');title.textContent=text;
  label.append(input,title);
  if (secondary) {const meta=document.createElement('span');meta.className='secondary';meta.textContent=secondary;label.append(meta);}
  container.append(label);
}
function drawCells() {
  $('cells').replaceChildren();
  if (!state.cells.length) {$('cells').textContent='No fort.6 files found.';return;}
  for (const cell of state.cells) {
    const label=document.createElement('label');label.className='cell-choice';
    const input=document.createElement('input');input.type='checkbox';input.value=cell.id;input.checked=state.chosen.has(cell.id);
    input.addEventListener('change',()=>{input.checked?state.chosen.add(cell.id):state.chosen.delete(cell.id);loadSelection();});
    label.append(input,document.createTextNode(cell.id));$('cells').append(label);
  }
}
async function discover(refresh=false) {
  const previous=new Set(state.chosen);
  const result=await api('/api/cells');state.cells=result.cells;
  $('project-directory').value=result.root;
  connectViewer(result.session_token);
  if (refresh) {state.cache.clear();state.chosen=new Set(state.cells.filter(c=>previous.has(c.id)).map(c=>c.id));}
  else if (state.cells.length) state.chosen.add(state.cells[0].id);
  $('local-status').textContent=`Local · ${state.cells.length} ${state.cells.length===1?'cell':'cells'}`;
  drawCells();await loadSelection();
}
async function loadSelection() {
  const request=++state.request;
  const missing=[...state.chosen].filter(id=>!state.cache.has(id));
  $('loading').hidden=!missing.length;
  if(missing.length) for(const id of ['csv','svg','png'])$(id).disabled=true;
  $('loading-text').textContent=`Reading ${missing.length} ${missing.length===1?'cell':'cells'}…`;
  const failures=[];
  // Two concurrent readers keep UI updates responsive with large result files.
  let next=0;
  async function worker() {
    while (next<missing.length) {
      const id=missing[next++];
      try {const cell=await api(`/api/cell/${encodeURIComponent(id)}`);if(request===state.request) state.cache.set(id,cell);}
      catch(e){failures.push(`${id}: ${e.message}`);}
    }
  }
  await Promise.all([worker(),worker()]);
  if (request!==state.request) return;
  $('loading').hidden=true;error(failures);
  updateFilters();updateChoices();render();
}
function fillFilter(id, values, current, allLabel) {
  const select=$(id);select.replaceChildren(new Option(allLabel,'all'));
  for (const value of [...values].sort((a,b)=>a-b)) select.add(new Option(String(value),String(value)));
  select.value=[...values].some(v=>String(v)===current)?current:'all';
  return select.value;
}
function updateFilters() {
  const samples=loaded().flatMap(c=>c.families[state.family]?.samples??[]);
  state.zone=fillFilter('zone',new Set(samples.map(s=>s.zone)),state.zone,'All');
  state.timeSet=fillFilter('time-set',new Set(samples.map(s=>s.set)),state.timeSet,'All');
}
function availableTimes() {
  const map=new Map();
  for(const cell of loaded()) for(const s of filteredSamples(cell.families.gamma,config())) map.set(String(s.time),s);
  return [...map.values()].sort((a,b)=>a.time-b.time);
}
function drawIsotopes() {
  const all=rankedSeries(loaded(),config());
  const names=new Set(all.map(s=>s.name));
  if(loaded().length) state.series=new Set([...state.series].filter(n=>names.has(n)));
  $('series-count').textContent=`${state.series.size} selected`;
  const search=$('isotope-search').value.toLowerCase();
  const choices=all.filter(s=>s.name.toLowerCase().includes(search)||friendly(s.name).toLowerCase().includes(search));
  // Selected items stay visible, then remaining matches ranked by peak value.
  choices.sort((a,b)=>Number(state.series.has(b.name))-Number(state.series.has(a.name)));
  $('isotopes').replaceChildren();
  for(const item of choices.slice(0,100)) option($('isotopes'),item.name,friendly(item.name),state.series.has(item.name),checked=>{
    checked?state.series.add(item.name):state.series.delete(item.name);drawIsotopes();render();
  },scientific(item.peak));
  if(choices.length>100){const p=document.createElement('p');p.className='hint';p.textContent=`${choices.length} matches. Search to find other isotopes.`;$('isotopes').append(p);}
}
function updateChoices() {
  drawAxisLimits();
  const gamma=state.family==='gamma';
  $('gamma-mode').hidden=!gamma;
  $('isotope-controls').hidden=gamma;
  $('gamma-times').hidden=!spectrum();
  $('gamma-bins').hidden=!gamma||spectrum();
  $('time-unit-wrap').hidden=spectrum();
  $('log-x').checked=axis().x;$('log-y').checked=axis().y;
  document.querySelectorAll('[data-mode]').forEach(b=>b.classList.toggle('active',b.dataset.mode===state.mode));
  const data=loaded().map(c=>c.families[state.family]).find(d=>d?.samples.length);
  $('unit-note').textContent=data?(gamma?'Volumetric source; no additional normalization.':state.family==='dose'?'Original ACAB surface dose.':'Results per unit volume.'):'No results loaded.';
  $('range-label').textContent=`Time filter [${units[state.timeUnit].label}]`;
  if (!gamma) {drawIsotopes();return;}
  const times=availableTimes();const values=new Set(times.map(s=>String(s.time)));
  state.times=new Set([...state.times].filter(v=>values.has(v)));
  if(!state.times.size&&times.length) state.times.add(String(times.at(-1).time));
  $('times').replaceChildren();
  for (const sample of times) option($('times'),String(sample.time),sample.label||'0 seconds',state.times.has(String(sample.time)),checked=>{
    checked?state.times.add(String(sample.time)):state.times.delete(String(sample.time));render();
  },`${scientific(sample.time)} s`);
  const bins=[...new Set(loaded().flatMap(c=>c.families.gamma.bins?.map(b=>b.upper)??[]))].sort((a,b)=>a-b);
  const validBins=new Set(['TOTAL',...bins.map(String)]);
  state.bins=new Set([...state.bins].filter(v=>validBins.has(v)));
  $('bins').replaceChildren();
  for (const value of ['TOTAL',...bins.map(String)]) option($('bins'),value,value==='TOTAL'?'Printed total':`≤ ${number(value)} MeV`,state.bins.has(value),checked=>{
    checked?state.bins.add(value):state.bins.delete(value);render();
  });
}
function empty(title,text) {
  $('curve-styles').hidden=true;
  $('chart').hidden=true;$('empty').hidden=false;$('empty-title').textContent=title;$('empty-text').textContent=text;
  $('table-view').hidden=true;
  for(const id of ['csv','svg','png']) $(id).disabled=true;
}
function updateSources(cells) {
  $('source-note').textContent=`${cells.map(c=>c.source).join(' · ')}${cells.length?' · Durations: 365-day year, 30-day month.':''}`;
  const warnings=cells.flatMap(c=>c.warnings.map(w=>`${c.cell}: ${w}`));
  $('warnings-wrap').hidden=!warnings.length;$('warnings-title').textContent=`Calculation warnings (${warnings.length})`;
  $('warnings').replaceChildren();for(const warning of warnings){const p=document.createElement('p');p.textContent=warning;$('warnings').append(p);}
}
function drawTable() {
  const rows=records(state.raw,state.family);
  $('table-summary').textContent=`${rows.length.toLocaleString('en-US')} observations · up to 500 shown · the CSV contains all, including zeros.`;
  const table=$('data-table');table.replaceChildren();if(!rows.length)return;
  const keys=Object.keys(rows[0]);const thead=document.createElement('thead');const tr=document.createElement('tr');
  for(const key of keys){const th=document.createElement('th');th.textContent=key;tr.append(th);}thead.append(tr);table.append(thead);
  const tbody=document.createElement('tbody');for(const row of rows.slice(0,500)){
    const tr=document.createElement('tr');for(const key of keys){const td=document.createElement('td');const v=row[key];td.textContent=typeof v==='number'&&!['zone','time_set'].includes(key)?scientific(v):String(v);if(typeof v==='number')td.className='numeric';tr.append(td);}tbody.append(tr);
  }table.append(tbody);
}
function render() {
  const cells=loaded();updateSources(cells);
  const spec=spectrum(),context=config();
  $('title').textContent=labels[state.family];
  const data=cells.map(c=>c.families[state.family]).find(d=>d?.samples.length);
  $('units').textContent=data?.unit??'—';
  $('axis-y-label').textContent=`Y [${data?.unit??'—'}]`;
  $('subtitle').textContent=spec?'Comparison of gamma sources by energy group.':'Time evolution of the selected series.';
  $('range-summary').textContent=`${cells.length} ${cells.length===1?'cell loaded':'cells loaded'}`;
  state.raw=makeTraces(cells,context);
  if(!cells.length)return empty(state.cells.length?'Select a cell':'Choose a data directory',state.cells.length?'Select one or more cells to load their results.':'Click Browse or enter a directory containing fort.6 or subdirectories with fort.6.');
  if(context.minimum!==null&&context.maximum!==null&&context.minimum>context.maximum)return empty('Invalid interval','The start time must be less than or equal to the end time.');
  if(state.raw.length>80){state.raw=[];return empty('Too many curves','Select fewer cells, series or times. The viewer supports up to 80 simultaneous curves.');}
  if(!state.raw.length)return empty('No series for this selection','Select an isotope, energy group or available time. Also check the zone and time set filters.');
  const prepared=maskForLog(state.raw,axis().x,axis().y);
  if(!prepared.visible)return empty('No visible points','The series have no positive values in the selected interval. Disable the logarithmic scale or change the filters.');
  const cellNames=[...new Set(prepared.series.map(t=>t.cell))];
  const seriesNames=[...new Set(prepared.series.map(t=>t.series))];
  $('curve-colors').replaceChildren();
  $('curve-styles').hidden=false;
  const plots=prepared.series.map(trace=>{
    const cellIndex=cellNames.indexOf(trace.cell),seriesIndex=seriesNames.indexOf(trace.series);
    const colorKey=JSON.stringify([state.family,spec,trace.name]);
    const color=state.curveColors.get(colorKey)??colors[(cellNames.length>1?cellIndex:seriesIndex)%colors.length];
    const label=document.createElement('label');label.className='curve-color';
    const picker=document.createElement('input');picker.type='color';picker.value=color;
    picker.setAttribute('aria-label',`Color for ${trace.name}`);
    picker.addEventListener('change',()=>{state.curveColors.set(colorKey,picker.value);render();});
    const text=document.createElement('span');text.textContent=trace.name;
    label.append(picker,text);$('curve-colors').append(label);
    const points=spec?spectrumSteps(trace.points,axis().x):trace.points;
    return {type:'scatter',mode:spec?'lines':'lines+markers',name:escape(trace.name),
      x:points.map(p=>p.x),y:points.map(p=>p.y),connectgaps:false,
      line:{color,width:1.9*(trace.series==='TOTAL'?1.5:1),dash:dashes[(cellNames.length>1?seriesIndex:0)%dashes.length],shape:'linear'},
      marker:{color,size:4},customdata:points.map(p=>[escape(p.label||'0 seconds'),p.time,p.zone,p.set,p.lower,p.upper]),
      hovertemplate:`<b>${escape(trace.name)}</b><br>Time: %{customdata[0]}<br>%{customdata[1]:.5e} s · Zone %{customdata[2]} · Set %{customdata[3]}<br>${spec?'Group: %{customdata[4]} – %{customdata[5]} MeV<br>':''}%{y:.5e} ${escape(trace.unit)}<extra></extra>`};
  });
  $('empty').hidden=true;$('chart').hidden=state.view!=='chart';$('table-view').hidden=state.view!=='table';
  for(const id of ['csv','svg','png'])$(id).disabled=false;
  $('subtitle').textContent=`${spec?'Spectrum by energy':'Time evolution'} · ${state.raw.length} ${state.raw.length===1?'curve':'curves'} · ${state.raw.reduce((n,t)=>n+t.points.length,0)} observations`;
  const notes=[];
  if(prepared.masked) notes.push(`${prepared.masked} non-positive points hidden by the logarithmic scale; preserved in the CSV.`);
  if(spec)notes.push('Rates integrated per group, without division by bin width.');
  notes.push('Drag to zoom; double-click to reset. Click the legend to hide a curve.');
  $('plot-note').textContent=notes.join(' ');
  const layout={margin:{l:87,r:28,t:28,b:100},paper_bgcolor:'#fff',plot_bgcolor:'#fff',
    font:{family:'Arial, sans-serif',size:state.fontSize,color:'#35545f'},hovermode:'closest',
    xaxis:{type:axis().x?'log':'linear',title:{text:spec?'Group upper boundary [MeV]':`Time [${units[state.timeUnit].label}]`,standoff:18},gridcolor:'#edf1f2',linecolor:'#cbd7db',ticks:'outside',zeroline:false,exponentformat:'power',automargin:true},
    yaxis:{type:axis().y?'log':'linear',title:{text:data.unit,standoff:15},gridcolor:'#e6edef',linecolor:'#cbd7db',zeroline:false,ticks:'outside',exponentformat:'power',automargin:true},
    legend:{orientation:'h',y:-.22,x:0,font:{size:state.legendSize}},
    uirevision:`${state.family}:${spec}:${state.timeUnit}:${state.zone}:${state.timeSet}:${state.minimum}:${state.maximum}:${JSON.stringify(limits())}:${axis().x}:${axis().y}`};
  const bounds=limits();
  Object.assign(layout.xaxis,axisRange(bounds.xmin===null?null:bounds.xmin/xFactor(),bounds.xmax===null?null:bounds.xmax/xFactor(),axis().x));
  Object.assign(layout.yaxis,axisRange(bounds.ymin,bounds.ymax,axis().y));
  Plotly.react($('chart'),plots,layout,{responsive:true,displaylogo:false,locale:'en',scrollZoom:true,
    modeBarButtonsToRemove:['lasso2d','select2d','sendChartToCloud'],toImageButtonOptions:{format:'png',filename:'apypa_results',scale:2}});
  if(state.view==='table')drawTable();
}

$('family').addEventListener('change',()=>{state.family=$('family').value;state.series=new Set(['TOTAL']);$('isotope-search').value='';updateFilters();updateChoices();render();});
$('isotope-search').addEventListener('input',drawIsotopes);
for(const [id,key] of [['font-size','fontSize'],['legend-size','legendSize']]) {
  const apply=()=>{if(!$(id).checkValidity()||$(id).value==='')return;state[key]=Number($(id).value);render();};
  $(id).addEventListener('input',apply);$(id).addEventListener('change',apply);
}
$('all-cells').onclick=()=>{state.chosen=new Set(state.cells.map(c=>c.id));drawCells();loadSelection();};
$('no-cells').onclick=()=>{state.chosen.clear();drawCells();loadSelection();};
$('only-total').onclick=()=>{state.series=new Set(['TOTAL']);drawIsotopes();render();};
$('no-series').onclick=()=>{state.series.clear();drawIsotopes();render();};
function applyTop() {
  const input=$('top-count'),count=Number(input.value);
  if(!input.reportValidity()||!Number.isSafeInteger(count)||count<1)return;
  const ranked=rankedSeries(loaded(),config()).filter(s=>s.name!=='TOTAL');
  state.series=new Set(['TOTAL',...ranked.slice(0,count).map(s=>s.name)]);
  drawIsotopes();render();
}
$('apply-top').onclick=applyTop;
$('top-count').addEventListener('keydown',event=>{if(event.key==='Enter')applyTop();});
for(const id of ['zone','time-set']) $(id).addEventListener('change',()=>{state[id==='zone'?'zone':'timeSet']=$(id).value;updateChoices();render();});
for(const id of ['minimum','maximum']) {
  const applyRange=event=>{
    if(!$(id).checkValidity()){if(event.type==='change')$(id).reportValidity();return;}
    state[id]=$(id).value===''?null:Number($(id).value)*units[state.timeUnit].factor;
    updateChoices();render();
  };
  $(id).addEventListener('input',applyRange);
  $(id).addEventListener('change',applyRange);
}
$('time-unit').addEventListener('change',()=>{state.timeUnit=$('time-unit').value;for(const id of ['minimum','maximum']) $(id).value=state[id]===null?'':state[id]/units[state.timeUnit].factor;updateChoices();render();});
function limitError(message) {$('axis-limit-error').textContent=message;$('axis-limit-error').hidden=false;}
function applyAxisLimits() {
  const values={};
  for(const key of ['xmin','xmax','ymin','ymax']) {
    if(!$(key).reportValidity())return;
    values[key]=$(key).value===''?null:Number($(key).value)*(key.startsWith('x')?xFactor():1);
  }
  try {axisRange(values.xmin,values.xmax,axis().x);axisRange(values.ymin,values.ymax,axis().y);}
  catch(e){limitError(e.message);return;}
  state.limits.set(limitKey(),values);$('axis-limit-error').hidden=true;render();
}
$('apply-axis-limits').onclick=applyAxisLimits;
for(const key of ['xmin','xmax','ymin','ymax'])$(key).addEventListener('keydown',event=>{if(event.key==='Enter')applyAxisLimits();});
for(const id of ['log-x','log-y'])$(id).addEventListener('change',()=>{
  const key=id==='log-x'?'x':'y',bounds=limits();
  try {axisRange(bounds[`${key}min`],bounds[`${key}max`],$(id).checked);}
  catch(e){$(id).checked=axis()[key];limitError(e.message);return;}
  axis()[key]=$(id).checked;$('axis-limit-error').hidden=true;render();
});
document.querySelectorAll('[data-mode]').forEach(button=>button.onclick=()=>{state.mode=button.dataset.mode;updateChoices();render();});
document.querySelectorAll('[data-view]').forEach(button=>button.onclick=()=>{state.view=button.dataset.view;document.querySelectorAll('[data-view]').forEach(b=>b.classList.toggle('active',b===button));render();});
$('last-time').onclick=()=>{const times=availableTimes();state.times=new Set(times.length?[String(times.at(-1).time)]:[]);updateChoices();render();};
$('three-times').onclick=()=>{const times=availableTimes();state.times=new Set(times.length?[times[0],times[Math.floor((times.length-1)/2)],times.at(-1)].map(s=>String(s.time)):[]);updateChoices();render();};
$('reset').onclick=()=>{state.limits.delete(limitKey());drawAxisLimits();render();Plotly.relayout($('chart'),{'xaxis.range':[null,null],'yaxis.range':[null,null],'xaxis.autorange':true,'yaxis.autorange':true});};
$('refresh').onclick=async()=>{try{await discover(true);}catch(e){error([e.message]);$('loading').hidden=true;}};
function download(blob,name){const a=document.createElement('a');const url=URL.createObjectURL(blob);a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
$('csv').onclick=()=>download(new Blob([csvText(records(state.raw,state.family))],{type:'text/csv;charset=utf-8'}),`apypa_${state.family}.csv`);
for(const format of ['svg','png'])$(format).onclick=()=>Plotly.downloadImage($('chart'),{format,filename:`apypa_${state.family}`,width:1500,height:850,scale:2});
window.addEventListener('resize',()=>{if($('chart').data)Plotly.Plots.resize($('chart'));});
try {await discover();}catch(e){error([e.message]);$('loading').hidden=true;empty('Could not load results','Check that the server is still running and click Refresh data.');}

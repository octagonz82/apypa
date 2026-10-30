// Pure data operations, shared by the browser and automated Node tests.
export const units = {
  seconds: {factor:1, label:'s'}, hours:{factor:3600,label:'hours'},
  days:{factor:86400,label:'days'}, years:{factor:31536000,label:'years'},
};
export const labels = {activity:'Activity',composition:'Composition',heat:'Thermal power',
  dose:'Surface gamma dose',gamma:'Gamma spectrum'};

export function axisRange(minimum,maximum,log=false) {
  for(const value of [minimum,maximum]) {
    if(value!==null&&(!Number.isFinite(value)||(log&&value<=0)))
      throw new Error(log?'Limits must be greater than zero on a logarithmic scale.':'Limits must be finite numbers.');
  }
  if(minimum!==null&&maximum!==null&&minimum>=maximum)
    throw new Error('The minimum must be less than the maximum.');
  const convert=value=>value===null?null:log?Math.log10(value):value;
  return {range:[convert(minimum),convert(maximum)],
    autorange:minimum===null?(maximum===null?true:'min'):(maximum===null?'max':false)};
}

export function filteredSamples(data, config) {
  return data.samples.map((s,index)=>({...s,index})).filter(s=>
    (config.zone==='all'||s.zone===Number(config.zone)) &&
    (config.timeSet==='all'||s.set===Number(config.timeSet)) &&
    (config.minimum===null||s.time>=config.minimum) &&
    (config.maximum===null||s.time<=config.maximum));
}

export function rankedSeries(cells, config) {
  const peaks=new Map();
  for (const cell of cells) {
    const data=cell.families[config.family];
    if (!data) continue;
    const indices=filteredSamples(data,config).map(s=>s.index);
    for (const series of data.series) {
      const peak=indices.reduce((v,i)=>Math.max(v,Math.abs(series.values[i]??0)),0);
      peaks.set(series.name,Math.max(peaks.get(series.name)??0,peak));
    }
  }
  return [...peaks].map(([name,peak])=>({name,peak})).sort((a,b)=>
    (a.name==='TOTAL'?-1:b.name==='TOTAL'?1:0)||b.peak-a.peak||a.name.localeCompare(b.name));
}

function evolutionGroups(samples) {
  const zones=new Map();
  for (const s of samples) {
    if (!zones.has(s.zone)) zones.set(s.zone,[]);
    zones.get(s.zone).push(s);
  }
  const groups=[];
  for (const [zone,points] of zones) {
    const split=new Set(points.map(s=>s.time)).size!==points.length;
    const sets=split?[...new Set(points.map(s=>s.set))]:[null];
    for (const set of sets) groups.push({zone,set,samples:points.filter(s=>set===null||s.set===set).sort((a,b)=>a.time-b.time)});
  }
  return groups;
}

export function makeTraces(cells, config) {
  const traces=[];
  const showZone=new Set(cells.flatMap(c=>c.families[config.family]?.samples.map(s=>s.zone)??[])).size>1;
  const factor=units[config.timeUnit].factor;
  for (const cell of cells) {
    const data=cell.families[config.family];
    if (!data) continue;
    const samples=filteredSamples(data,config);
    if (config.family==='gamma'&&config.mode==='spectrum') {
      for (const sample of samples.filter(s=>config.times.has(String(s.time)))) {
        const points=data.bins.map((bin,i)=>({x:bin.upper,y:data.values[sample.index][i],
          ...bin,time:sample.time,label:sample.label,set:sample.set,zone:sample.zone}));
        traces.push({cell:cell.cell,series:sample.label,unit:data.unit,
          name:`${cell.cell} · ${sample.label}${showZone?` · Z${sample.zone}`:''} · S${sample.set}`,points,spectrum:true});
      }
    } else {
      const selected=config.family==='gamma'?[...config.bins]:[...config.series];
      for (const name of selected) {
        let series=data.series.find(s=>s.name===name);
        let bin=null;
        if (config.family==='gamma') {
          const index=data.bins.findIndex(b=>String(b.upper)===name);
          if (name==='TOTAL') series={name,values:data.totals};
          else if (index!==-1) {
            bin=data.bins[index];
            series={name:`E ≤ ${bin.upper} MeV`,values:data.values.map(row=>row[index])};
          }
        }
        if (!series) continue;
        for (const group of evolutionGroups(samples)) {
          const points=group.samples.map(s=>({x:s.time/factor,y:series.values[s.index],
            time:s.time,label:s.label,set:s.set,zone:s.zone,lower:bin?.lower??null,upper:bin?.upper??null}));
          if (!points.length) continue;
          traces.push({cell:cell.cell,series:series.name,unit:data.unit,
            name:`${cell.cell} · ${series.name}${showZone?` · Z${group.zone}`:''}${group.set!==null?` · S${group.set}`:''}`,
            points,spectrum:false});
        }
      }
    }
  }
  return traces;
}

export function maskForLog(traces, logX, logY) {
  let masked=0,visible=0;
  const series=traces.map(trace=>({...trace,points:trace.points.map(point=>{
    if (point.y===null||point.y===undefined) return {...point,y:null};
    if ((logX&&point.x<=0)||(logY&&point.y<=0)) {masked++;return {...point,y:null};}
    visible++;return {...point};
  })}));
  return {series,masked,visible};
}

export function spectrumSteps(points, logX) {
  // Each bin has its own two edges: a zero/absent neighbour must not erase the
  // horizontal segment of a positive bin on a logarithmic Y axis.
  return points.flatMap(p=>(!logX||p.lower>0)?[{...p,x:p.lower},{...p}]:[{...p}]);
}

export function records(traces, family) {
  return traces.flatMap(trace=>trace.points.map(point=>({
    cell:trace.cell,zone:point.zone,time_set:point.set,
    time_s:point.time,time_text:point.label,quantity:labels[family],
    series:trace.series,unit:trace.unit,lower_energy_MeV:point.lower??'',
    upper_energy_MeV:point.upper??'',value:point.y??'',
  })));
}

export function csvText(rows) {
  if (!rows.length) return '';
  const escape=value=>{
    const text=String(value??'');
    // Keep identifiers that resemble Excel formulas literal in CSV imports.
    const literal=typeof value==='string'&&/^[=+@-]/.test(text)?`'${text}`:text;
    return /["\r\n,]/.test(literal)?`"${literal.replaceAll('"','""')}"`:literal;
  };
  const keys=Object.keys(rows[0]);
  return '\uFEFF'+[keys.map(escape).join(','),...rows.map(row=>keys.map(k=>escape(row[k])).join(','))].join('\r\n');
}

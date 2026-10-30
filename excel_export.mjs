// Excel rendering/export engine for acab_reader.py. No transformations of source values.
import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const input = JSON.parse(await fs.readFile(process.argv[2], 'utf8'));
const workbook = Workbook.create();
const previews = [];
const showZone = input.zones.length > 1;
for (const spec of input.sheets) {
  const sheet = workbook.worksheets.add(spec.name);
  sheet.showGridLines = false;
  sheet.getRange('A1').values = [[`Cell ${input.cell} · ${spec.name}`]];
  sheet.getRange('A1').format.font = {name: 'Arial', size: 14, bold: true};
  sheet.getRange('A2').values = [[`Source: ${spec.source}`]];
  sheet.getRange('A3').values = [[`Units: ${spec.unit}${spec.family === 'gamma' ? ` · ${spec.energy_header}` : ''}`]];
  const volume = (spec.volumes ?? []).filter(v => v !== null).join(', ');
  sheet.getRange('A4').values = [[
    spec.family === 'dose' ? `${(spec.dose_descriptions ?? []).join(' / ')} · Original dose in Sv/h.` :
    spec.family === 'gamma' ? 'Volumetric source per group; original ACAB values.' :
    `Isotopes retained: ${spec.retained} · Omitted as all zero: ${spec.omitted}`
  ]];
  const divided = ['composition', 'activity', 'heat'].includes(spec.family);
  sheet.getRange('A5').values = [[`Volume [cm³]: ${volume || 'not printed'} · ${divided ? 'Results / volume' : 'No additional volume division'} · Year: 365 days; month: 30 days.`]];
  sheet.getRange('A2:A5').format.font = {name: 'Arial', size: 10, color: '#526174'};
  // Keep the reader's full schema; omit redundant zone labels only in Excel.
  const matrix = showZone ? spec.matrix : spec.family === 'gamma'
    ? spec.matrix.map(row => row.filter((_, i) => i !== 1))
    : spec.matrix.filter((_, i) => i !== 1);
  if (!matrix.length) {
    sheet.getRange('A7').values = [['No results with numeric times in this family.']];
    continue;
  }
  const n = matrix.length, m = matrix[0].length;
  if (n + 5 > 1048576 || m > 16384) throw new Error(`Worksheet ${spec.name} exceeds Excel limits.`);
  if (matrix.some(row => row.length !== m)) throw new Error(`Irregular matrix: ${spec.name}`);
  const used = sheet.getRangeByIndexes(5, 0, n, m);
  used.values = matrix;
  used.format.font = {name: 'Arial', size: 10};
  used.format.rowHeight = 19;
  used.format.columnWidth = 16;
  used.format.verticalAlignment = 'center';
  sheet.getRangeByIndexes(5, 0, n, 1).format.columnWidth = spec.family === 'gamma' ? 20 : 28;
  const headerRows = spec.family === 'gamma' ? 1 : showZone ? 4 : 3;
  const headers = sheet.getRangeByIndexes(5, 0, headerRows, m);
  headers.format.fill = '#233B55';
  headers.format.font = {name: 'Arial', size: 10, color: '#FFFFFF', bold: true};
  headers.format.horizontalAlignment = 'center';
  headers.format.borders = {insideVertical: {style: 'thin', color: '#FFFFFF'}};
  headers.format.rowHeight = spec.family === 'gamma' ? 34 : 24;
  headers.format.wrapText = true;
  if (spec.family === 'gamma') {
    const timeColumn = showZone ? 2 : 1, textColumn = timeColumn + 1, energyColumn = textColumn + 1;
    if (showZone) sheet.getRangeByIndexes(5, 1, n, 1).format.columnWidth = 9;
    sheet.getRangeByIndexes(5, timeColumn, n, 1).format.columnWidth = 26;
    sheet.getRangeByIndexes(5, textColumn, n, 1).format.columnWidth = 55;
    sheet.getRangeByIndexes(6, textColumn, n - 1, 1).format.wrapText = true;
    sheet.getRangeByIndexes(6, 0, n - 1, m).format.rowHeight = 32;
    sheet.getRangeByIndexes(5, m - 2, n, 2).format.columnWidth = 32;
    sheet.getRangeByIndexes(5, energyColumn, 1, m - energyColumn - 2).setNumberFormat('0.000E+00');
    sheet.getRangeByIndexes(6, timeColumn, n - 1, 1).setNumberFormat('0.000E+00');
    sheet.getRangeByIndexes(6, energyColumn, n - 1, m - energyColumn).setNumberFormat('0.000E+00');
  } else {
    sheet.getRangeByIndexes(5 + headerRows - 2, 1, 1, m - 1).setNumberFormat('0.000E+00');
    sheet.getRangeByIndexes(5 + headerRows, 1, n - headerRows, m - 1).setNumberFormat('0.000E+00');
    sheet.getRangeByIndexes(5 + headerRows - 1, 0, 1, m).format.rowHeight = 72;
    for (let r = headerRows; r < n; r++) {
      if (matrix[r][0] === 'TOTAL') {
        const total = sheet.getRangeByIndexes(r + 5, 0, 1, m);
        total.format.fill = '#E5ECF3';
        total.format.font.bold = true;
      }
    }
  }
  sheet.freezePanes.freezeRows(5 + headerRows);
  sheet.freezePanes.freezeColumns(spec.family === 'gamma' ? (showZone ? 3 : 2) : 1);
  previews.push({name: spec.name, range: 'A1:H14'});
}
workbook.recalculate();
const scan = await workbook.inspect({kind:'match', searchTerm:'#REF!|#DIV/0!|#VALUE!|#NUM!',
  options:{useRegex:true, maxResults:10}, maxChars:1000});
if (input.preview_dir) {
  await fs.mkdir(input.preview_dir, {recursive:true});
  for (const item of previews) {
    const preview = await workbook.render({sheetName:item.name, range:item.range, scale:1.5, format:'png'});
    await fs.writeFile(path.join(input.preview_dir, `${item.name}.png`), new Uint8Array(await preview.arrayBuffer()));
  }
  const inspection = await workbook.inspect({kind:'region', sheetId:'Activity', range:'A6:E12', maxChars:1800});
  await fs.writeFile(path.join(input.preview_dir, 'inspection.json'), JSON.stringify({scan:scan.ndjson, sample:inspection.ndjson}));
}
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(input.output);

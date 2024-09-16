updateViz = function (id) {
  index = localStorage.getItem("genome") ? localStorage.getItem("genome") : 0;
  indexSelected = index;
  genome_selected = genomes[index];
  $("#vizStructure").attr("src", `../SecondaryStructure/${genome_selected}/Visualizations/${id}`);
  $("#collapseExample").collapse("show");
  $("html,body").animate(
    {
      scrollTop: $("#collapseExample").offset().top,
    },
    "slow"
  );
};

let sum = function (array) {
  return array.reduce(function (pv, cv) {
    return pv + cv;
  }, 0);
};
updateSummary = function (genome) {
  let valuesSummary = [];
  let data = shortDescription[genome];
  let naValues = [];
  for (let i of ["antisense", "orphan", "internal"]) {
    let subdata = data[i] || {};
    for (let j of ["RNA", "COD"]) {
      valuesSummary.push(subdata[j] || 0);
    }

    naValues.push(subdata["NA"] || 0);
  }
  valuesSummary.push(data["Ignored"]);
  valuesSummary.push(sum(valuesSummary));
  valuesSummary.push(sum(valuesSummary.slice(0, 2)));
  valuesSummary.push(sum(valuesSummary.slice(2, 4)));
  valuesSummary.push(sum(valuesSummary.slice(4, 6)));
  let ids = ["ncATSS", "CodingATSS", "ncOTSS", "CodingOTSS", "ncINTTSS", "CodingINTTSS", "iTSS", "allTSS", "ATSS", "OTSS", "INTTSS"];
  valuesSummary.forEach(function (value, i) {
    $(`#${ids[i]}`).text(value);
  });
  if (sum(naValues) > 0) {
    $("#noclassAS").text(naValues[0]);
    $("#ATSS").text(valuesSummary[8] + naValues[0]);
    $("#noclassOR").text(naValues[1]);
    $("#OTSS").text(valuesSummary[9] + naValues[1]);
    $("#noclassINT").text(naValues[2]);
    $("#INTTSS").text(valuesSummary[10] + naValues[2]);
    $("#allTSS").text(valuesSummary[7] + sum(naValues));
  } else {
    $("#not-classified-or").remove();
    $("#not-classified-as").remove();
    $("#not-classified-int").remove();
  }
};
createOption = function (text, value) {
  return $("<option>").val(value).text(text);
};

genomes.forEach(function (x, i) {
  $("#selectGenome").append(createOption(x, i));
});

let index = localStorage.getItem("genome") ? localStorage.getItem("genome") : 0;
indexSelected = index;
genome_selected = genomes[index];

$(`#selectGenome option[value='${indexSelected}']`).attr("selected", true);
$(`.link-to-motif`).each(function () {
  $(this).attr("href", `../MotifAnalysis/${genome_selected}/meme_out/meme.html`);
});
$("#selectGenome").on("change", function () {
  let index = $(`#selectGenome option`).filter(":selected").val();
  genome_selected = genomes[index];
  console.log("genome_selected:", genome_selected);
  localStorage.setItem("genome", index);
  let currentTable = $(".dataTables_wrapper")[0].id.replace("_wrapper", "");
  let data;
  switch (currentTable) {
    case "dataTableOverview":
      data = overviewData[genome_selected];
      break;
    case "dataTableClass":
      data = classified[genome_selected];
      break;
    case "dataTableTerms":
      data = terminators[genome_selected];
      break;
    case "dataTableAvoided":
      data = avoidedTSS[genome_selected];
      break;
  }
  let datatable = new $.fn.dataTable.Api(`#${currentTable}`);
  datatable.clear();
  datatable.rows.add(data);
  if (currentTable == "dataTableOverview") {
    dataSummary = summaryMotifs[genome_selected];
    newSummary = dataSummary
      .map((x) => `${x[0]} (${x[1]})<br/>Appears ${x[3]} times (E-value: ${x[2].toExponential(2)})`)
      .join("<br/>");
    $(datatable.column(9).header()).html(
      '<a target="_blank" rel="noopener noreferrer" data-html="true" class="wide-tooltip link-to-motif" data-toggle="tooltip" title="' +
        newSummary +
        '">Promoter motifs</a>'
    );
  }
  datatable.draw();

  $(`.link-to-motif`).each(function () {
    $(this).attr("href", `../MotifAnalysis/${genome_selected}/meme_out/meme.html`);
  });
  updateSummary(genome_selected);
});

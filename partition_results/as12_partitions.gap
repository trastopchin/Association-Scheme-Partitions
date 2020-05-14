for i in [2..Length(as12)] do
  LogPartitionsScheme("./partition_results/as12_partitions.gap", as12[i], Concatenation("as12_", String(i)));
od;

Read("equitable_partitions.gap");

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("./hanaki_association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes of order 5 and 16
# as05 is a list of the schemes of order 5
Read("./classified_schemes/schemes_order05.gap");
Read("./classified_schemes/schemes_order06.gap");
Read("./classified_schemes/schemes_order07.gap");
Read("./classified_schemes/schemes_order08.gap");
Read("./classified_schemes/schemes_order09.gap");
Read("./classified_schemes/schemes_order10.gap");
Read("./classified_schemes/schemes_order11.gap");
Read("./classified_schemes/schemes_order12.gap");
#
# Print("-----------------------------------\n");
# Print("Schemes of order 5:\n");
# for i in [2..Length(as05)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as05[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 6:\n");
# for i in [2..Length(as06)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as06[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 7:\n");
# for i in [2..Length(as07)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as07[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 8:\n");
# for i in [2..Length(as08)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as08[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 9:\n");
# for i in [2..Length(as09)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as09[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 10:\n");
# for i in [2..Length(as10)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as10[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 11:\n");
# for i in [2..Length(as11)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as11[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 12:\n");
# for i in [2..Length(as12)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   EquitablePartitions(as12[i]);
# od;




# for i in [2..Length(as05)] do
#   LogPartitionsScheme("./partition_results/as05_partitions", as05[i], Concatenation("as05_", String(i)));
# od;
#
# for i in [2..Length(as06)] do
#   LogPartitionsScheme("./partition_results/as06_partitions", as06[i], Concatenation("as06_", String(i)));
# od;
#
# for i in [2..Length(as07)] do
#   LogPartitionsScheme("./partition_results/as07_partitions", as07[i], Concatenation("as07_", String(i)));
# od;
#
# for i in [2..Length(as08)] do
#   LogPartitionsScheme("./partition_results/as08_partitions", as08[i], Concatenation("as08_", String(i)));
# od;
#
# for i in [2..Length(as09)] do
#   LogPartitionsScheme("./partition_results/as09_partitions", as09[i], Concatenation("as09_", String(i)));
# od;
#
# for i in [2..Length(as10)] do
#   LogPartitionsScheme("./partition_results/as10_partitions", as10[i], Concatenation("as10_", String(i)));
# od;
#
# for i in [2..Length(as11)] do
#   LogPartitionsScheme("./partition_results/as11_partitions", as11[i], Concatenation("as11_", String(i)));
# od;

for i in [2..Length(as12)] do
  LogPartitionsScheme("./partition_results/as12_partitions", as12[i], Concatenation("as12_", String(i)));
od;

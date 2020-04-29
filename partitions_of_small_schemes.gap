Read("enumerate_partitions.gap");

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("./association_scheme.gap");

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

Print("-----------------------------------\n");
Print("Schemes of order 5:\n");
for i in [2..Length(as05)] do
  Print("\nScheme No.");
  Print(i);
  Print(":");
  ComputeGoodPartitions(as05[i]);
od;

Print("-----------------------------------\n");
Print("Schemes of order 6:\n");
for i in [2..Length(as06)] do
  Print("\nScheme No.");
  Print(i);
  Print(":");
  ComputeGoodPartitions(as06[i]);
od;

# Print("-----------------------------------\n");
# Print("Schemes of order 7:\n");
# for i in [2..Length(as07)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   ComputeGoodPartitions(as07[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 8:\n");
# for i in [2..Length(as08)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   ComputeGoodPartitions(as08[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 9:\n");
# for i in [2..Length(as09)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   ComputeGoodPartitions(as09[i]);
# od;
#
# Print("-----------------------------------\n");
# Print("Schemes of order 10:\n");
# for i in [2..Length(as10)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   ComputeGoodPartitions(as10[i]);
# od;

# Print("-----------------------------------\n");
# Print("Schemes of order 11:\n");
# for i in [2..Length(as11)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   ComputeGoodPartitions(as11[i]);
# od;

# Print("-----------------------------------\n");
# Print("Schemes of order 12:\n");
# for i in [2..Length(as12)] do
#   Print("\nScheme No.");
#   Print(i);
#   Print(":");
#   ComputeGoodPartitions(as12[i]);
# od;
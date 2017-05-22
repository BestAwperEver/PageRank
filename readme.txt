Command line parameters:
filename=X precision=X [test=[true|false]] [dumping_factor=X]

Default values for test and for dumping_factor are false and .85, respectivly

Test file format:
Number_of_sites
Donor_site_name Recipient1 Recipient2 ...
...

Example:
8
Home About Products Links
About Home
Products Home
Links Home ExtA ExtB ExtC ExtD
ExtA Home
ExtB Home
ExtC Home
ExtD Home

"Real data" file format:
Donor_index1 Recipient_index1
Donor_index2 Recipient_index2
...

Example:
1	2
1	3
1 	4
2	4
2	5
...

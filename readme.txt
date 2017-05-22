Command line parameters:
filename=X				// file with the data
precision=X				// sensitivity of the algorithm
[output_filename=X]			// file with output weights; isn't created, if not specified
[test={true|false}]			// whether data is in test format; default false
[dumping_factor=X]			// dumping factor for the algorithm; default 0.85
[power_iteration={true|false}]		// use base Power Iteration algorithm for obtaining a baseline; default false

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

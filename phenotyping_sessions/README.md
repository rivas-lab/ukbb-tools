# Phenotyping Sessions

This guide will take you through how best to define phenotypes as a group.

## What is a phenotyping session? Why do we do it?

A phenotyping session is what we should convene when new data is added to our existing dataset. Decisions made at a phenotyping session allow the lab to extract binary or quantitative (QT) phenotypes from new or updated BioBank fields, and subsequently allow the lab to run a number of downstream analyses including GWAS and PheWAS.

## Getting the group together

This might just be the hardest part of this section - you'll need to find a time when all of the group is ready, willing, and able to work on this together. Doing it by your lonesome is sucky and should never happen. Ping the Slack well in advance and get a quorum of people together before these things happen.

## At the session

At the beginning of the session, there will be a `.tsv` on Google Drive shared with the group, formatted like so:

BLAH

The left-most column is the Field ID from BioBank, _____________ETC_____________.

At the BioBank phenotyping session, the concept is to divide and conquer. You'll be given a number of phenotypes (field IDs) to define. Look up the field ID [in the BioBank search box](http://biobank.ctsu.ox.ac.uk/crystal/search.cgi). This should give you the name of the field ID when you click on the corresponding number. Fill this name in in the appropriate column in the Google sheet. Then, for each row, you should ask yourself:

1) Is this field referring to one trait or multiple? 
- One trait would be something like "Have you smoked before in your life?"
- Multiple traits would be something like "Please mark if you have had stroke, cardiac event, both, or none" (you would be splitting this into two phenotypes - "Had stroke" and "Had cardiac event").
    - In the case of multiple-trait fields, you will need to copy paste that row that you're on as many times as needed (in the case above, copy and paste it once to have a total of 2 rows). 
  
2) Assuming that the row now refers to a single trait and not multiple, is this a binary trait or a quantitative (QT) trait? 
- In the case of it being a binary trait: your decisions are pretty simple. Find out what value corresponds to cases and what value corresponds to controls via the [BioBank search box](http://biobank.ctsu.ox.ac.uk/crystal/search.cgi).
  - Specifically, after searching for the field ID, when you click on the field ID in question, you should see some tabs, with "Data" being highlighted by default. Right under that, there should be some text saying "X items of data are available, covering Y participants, encoded using Data-Coding Z.") Click on the number "Z" to find out how the data is encoded.
  - According to the convention that we have, anything that is not case or control will automatically be encoded as missing.
- In the case of it being a quantitative trait, we need to look at how it is encoded. Some quantitative traits are encoded just as is, like [age started smoking in current smokers](http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=3436), which is encoded in the BioBank as years. If this is the case, leave all fields except name and field ID blank, but also include the special encodings into the Google sheet (for the link above, -1 means "Do not know" and -3 means "Prefer not to answer").
  - Some quantitative traits were measured using multiple choice questions. Imagine that you have a question like the following:
    
    Which of the following quantities typifies the quantity of sleep you get?
    - 0-2 hours
    - 2-4 hours
    - 4-6 hours
    - 6-8 hours
    - 8-10 hours
    - 10+ hours
    
    For this question, you should split your existing row into 6 rows. The naming convention for these kinds of rows is "Field Description (Parenthetical additional descriptor)" (e.g. "Average amount of sleep (0-2 hours)"). Additionally, there is an "ORDER" column in the sheet. For this kind of quantitative trait, the person doing this trait should ordinally assign values to this column based on quantity (1 for "0-2 hours", 2 for "2-4 hours", etc.).
  
  # TODO: For QT traits, if there are special values in the field that we need to remap to other ones (e.g. for some dietary traits, -10 means less than 1), how do we do this?
  
## After the session
  
You are ready to now compile `.phe` files and [run QC](https://github.com/rivas-lab/ukbb-tools/tree/master/qc).

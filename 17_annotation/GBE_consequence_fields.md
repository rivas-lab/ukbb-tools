# GBE consequence fields

When loading the variant annotation to GBE, we convert the variant annotation string so that it is loadble by GBE.

The schema of the annotation string is defined in the parser script `parsing.py` in `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser`.

(you can search it with `grep Allele *.py` in that directory)

We extracted the list of columns in [`GBE_consequence_fields.txt`](GBE_consequence_fields.txt).


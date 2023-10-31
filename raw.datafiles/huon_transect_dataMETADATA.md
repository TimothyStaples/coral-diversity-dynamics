--Metadata for huon_transect_data.csv--

This .csv provides raw transect data from uplifted Holocene coral reef terraces along the Huon Peninsula, Papua New Guinea. Subsets of these data were used to generate analyses, with full processing and analysis pathway in the top-level R-script.

Details of methodology are available at the publication connected to this repository, with additional details at:

Pandolfi, J. M., Tudhope, A. W., Burr, G., Chappell, J., Edinger, E., Frey, M., Steneck, R., Sharma, C., Yeates, A., Jennions, M., Lescinsky, H., & Newton, A. (2006). Mass mortality following disturbance in Holocene coral reefs from Papua New Guinea. Geology, 34(11), 949–952. https://doi.org/10.1130/G22814A.1

Edinger, E. N., Burr, G. S., Pandolfi, J. M., & Ortiz, J. C. (2007). Age accuracy and resolution of Quaternary corals used as proxies for sea level. Earth and Planetary Science Letters, 253(1), 37–49. https://doi.org/https://doi.org/10.1016/j.epsl.2006.10.014

Column metadata:

transect - (category) ID number of horizontal reef terrace transect.
country - (category) country location of transect (all = Papua New Guinea).
region - (category) region location of transect (all = Huon peninsula).
age - (category) epoch of transect (all = Holocene).
site - (category) spatial cluster of transect (referred to as "locality" in publication associated with this repository).
locality - (category) sample site of transect (referred to as "site" in publication associated with this repository). All transects at a given site are horizontally aligned, and separated only vertically
date - (date) date of sample
latitude - (numeric, digital degrees) latitude of transect (all transects at a given site share the same coordinates)
longitude - (numeric, digital degrees) longitude of transect (all transects at a given site share the same coordinates)
elevation - (free text) description of transect elevation relative to other features and transects.
height.from.dist - (numeric, metres) height of transect from a known disturbance layer.
height.ord - (category) vertical rank of transects up each site (1 = lowest, oldest transect).
dist.from.top - (numeric, metres) height of transect down from the top of the Holocene reef terrace.
age.min - (numeric, years before 1950AD) minimum calibrated transect age from 14C radiometric dates of coral fragments. See associated publication for chronology model methodology.
age.max - (numeric, years before 1950AD) maximum calibrated transect age from 14C radiometric dates of coral fragments. See associated publication for chronology model methodology.
age.median - (numeric, years before 1950AD) median calibrated transect age from 14C radiometric dates of coral fragments. See associated publication for chronology model methodology.
age.mean - (numeric, years before 1950AD) mean calibrated transect age from 14C radiometric dates of coral fragments. See associated publication for chronology model methodology.
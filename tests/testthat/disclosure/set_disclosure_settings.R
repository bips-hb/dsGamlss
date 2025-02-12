#
# DataSHIELD disclosure settings
#

set.standard.disclosure.settings <- function() {
    options(datashield.privacyLevel = "5")
    options(default.datashield.privacyControlLevel = "permissive")
    options(default.nfilter.glm = "0.33")
    options(default.nfilter.kNN = "3")
    options(default.nfilter.string = "80")
    options(default.nfilter.subset = "3")
    options(default.nfilter.stringShort = "20")
    options(default.nfilter.tab = "3")
    options(default.nfilter.noise = "0.25")
    options(default.nfilter.levels.density = "0.33")
    options(default.nfilter.levels.max = "40")
}

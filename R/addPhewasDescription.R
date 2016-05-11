addPhewasDescription <- function(data, keep.unmatched.rows=F,for.plots=F) {
  .Deprecated("addPhecodeInfo",msg="addPhewasDescription has been deprecated and will call addPhecodeInfo.
              Please see ?addPhecodeInfo for a description of methods applied.")
  addPhecodeInfo(data,groups=F,groupnums=F)
}
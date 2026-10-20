.anonymous_sample_registry <- new.env(parent = emptyenv())
anonymous_cohorts <- function(x) {
  values <- as.character(x)
  if (anyNA(values) || any(!grepl("^(?:Cohort_[A-Z]+(?:_(?:DLBCL|FL|MCL|BL))?|Training|Test|All)$", values, perl = TRUE)))
    stop("Cohort labels must use neutral Cohort aliases, Training, Test, or All.", call. = FALSE)
  values
}
anonymous_ids <- function(x, field = "identifier") {
  values <- as.character(x)
  if (anyNA(values) || any(!grepl("^[1-9][0-9]*$", values)))
    stop(paste0(field, " must contain pre-anonymized positive integer identifiers."), call. = FALSE)
  values
}
validate_inclusion <- function(x) {
  if (is.logical(x) && !anyNA(x)) return(x)
  values <- tolower(as.character(x))
  if (anyNA(values) || any(!values %in% c("true", "false", "1", "0")))
    stop("include_in_analysis must contain explicit boolean values.", call. = FALSE)
  values %in% c("true", "1")
}
validate_anonymous_data <- function(data) {
  if (!is.data.frame(data)) stop("Expected a data frame.", call. = FALSE)
  normalized <- tolower(gsub("[^a-zA-Z0-9]", "", names(data)))
  private_fields <- c("patientname", "subjectname", "participantname", "donorname", "fullname", "firstname", "lastname", "givenname", "familyname", "surname", "initials", "email", "emailaddress", "phone", "phonenumber", "telephone", "mobile", "address", "homeaddress", "streetaddress", "postaladdress", "postcode", "zipcode", "dateofbirth", "birthdate", "dob", "nationalid", "idcard", "passport", "ssn", "medicalrecordnumber", "medicalrecordid", "mrn", "hospitalnumber", "hospitalid", "originalpatientid", "originalsubjectid", "originalsampleid", "serumsampleid", "urinesampleid", "serumnumber", "urinenumber", "localpath", "filepath", "sourcepath", "rawpath", "notes", "note", "remarks", "remark", "comments", "comment")
  if (any(normalized %in% private_fields))
    stop("Input contains prohibited identifier or free-text fields.", call. = FALSE)
  if (any(grepl("\\p{Han}", names(data), perl = TRUE)))
    stop("Input column names must use the English analysis schema.", call. = FALSE)
  identifier_fields <- c("sample", "sampleid", "patient", "patientid", "subject", "subjectid", "donor", "donorid", "specimen", "specimenid", "patientkey", "cellid", "cellbarcode", "barcode", "spotid", "roiid", "aoiid", "fovid")
  for (field in names(data)[normalized %in% identifier_fields])
    anonymous_ids(data[[field]], field)
  for (field in names(data)[normalized %in% c("dataset", "datasetid", "cohort", "cohortid", "study", "studyid", "datasetlabel", "cohortlabel", "studylabel")])
    anonymous_cohorts(data[[field]])
  if ("samples" %in% names(data)) {
    values <- as.character(data[["samples"]])
    if (anyNA(values) || any(!grepl("^[1-9][0-9]*(\\+[1-9][0-9]*)*$", values)))
      stop("samples must contain numeric specimen identifiers separated by plus signs.", call. = FALSE)
  }
  for (field in names(data)) {
    x <- data[[field]]
    if (!is.character(x) && !is.factor(x)) next
    values <- as.character(x[!is.na(x)])
    if (any(grepl("\\p{Han}", values, perl = TRUE)))
      stop("Input character fields must use English analysis values.", call. = FALSE)
    if (any(grepl("^(?:/|~[/\\\\]|[A-Za-z]:[/\\\\]|\\\\\\\\)", values, perl = TRUE)))
      stop("Input contains local filesystem locations.", call. = FALSE)
  }
  if (all(c("sample_id", "patient_id") %in% names(data))) {
    samples <- anonymous_ids(data[["sample_id"]], "sample_id")
    patients <- anonymous_ids(data[["patient_id"]], "patient_id")
    links <- unique(data.frame(sample_id = samples, patient_id = patients))
    if (anyDuplicated(links$sample_id))
      stop("A sample_id is linked to more than one patient_id.", call. = FALSE)
    for (i in seq_len(nrow(links))) {
      key <- links$sample_id[i]
      patient <- links$patient_id[i]
      if (exists(key, envir = .anonymous_sample_registry, inherits = FALSE) &&
          !identical(get(key, envir = .anonymous_sample_registry, inherits = FALSE), patient))
        stop("Cross-file sample-to-patient linkage is inconsistent.", call. = FALSE)
      assign(key, patient, envir = .anonymous_sample_registry)
    }
  }
  data
}
source_release_dependency <- function(path) {
  if (!file.exists(path))
    stop(paste0("Required shared code is absent from the supplied archive: ", basename(path)), call. = FALSE)
  sys.source(path, envir = parent.frame())
  invisible(NULL)
}

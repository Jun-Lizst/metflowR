# checkData(data = grep("batch", x = temp.file, value = TRUE),
#           sample.info = "sample.info.csv",
#           path = ".")


setGeneric(name = "checkData",
           def = function(data = c("batch1.demo.csv", "batch2.demo.csv"),
                          sample.info = "sample.info.demo.csv",
                          path = "."){

             sample.info.record <- NULL
             cat("Read data.\n")
             cat("--------------------------------------------------------------\n")
             data <- lapply(data, function(x){
               as.data.frame(readr::read_csv(file.path(path, x),
                                             col_types = readr::cols(),
                                             progress = TRUE))
             })

             sample.info <- readr::read_csv(file.path(path, sample.info), col_types = readr::cols(),
                                            progress = TRUE)

             sample.info <- as.data.frame(sample.info)
             #####

             ##check data, component
             data.record <- lapply(seq_along(data), function(x){
               temp.data.record <- NULL
               data.col.name <- colnames(data[[x]])
               if(all(data.col.name != "name")){
                 cat("Error: No name in data.\n")
                 temp.data.record <- c(temp.data.record, "Error")
               }else{
                 # cat("OK: The first column of data is name.\n")
                 temp.data.record <- c(temp.data.record, "OK")
               }

               if(all(data.col.name != "mz")){
                 cat("Error: No mz in data.\n")
                 temp.data.record <- c(temp.data.record, "Error")
               }else{
                 # cat("OK: The second column of data is mz.\n")
                 temp.data.record <- c(temp.data.record, "OK")
               }

               if(all(data.col.name != "rt")){
                 cat("Error: No rt in data.\n")
                 temp.data.record <- c(temp.data.record, "Error")
               }else{
                 # cat("OK: The third column of data is not rt.\n")
                 temp.data.record <- c(temp.data.record, "OK")
               }
               temp.data.record
             })

             cat("--------------------------------------------------------------\n")
             ##check sample.info
             if(sum(is.na(sample.info)) > 0){
               cat("Error: There are", sum(is.na(sample.info)), "NAs in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               # cat("OK: There are no NAs in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "OK")
             }

             if(ifelse(is.na(sum(sample.info == "") > 0), FALSE, sum(sample.info == "") > 0)){
               cat("Error: There are", sum(sample.info == ""), "spaces in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               # cat("OK: There are no spaces in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "OK")
             }


             if(colnames(sample.info)[1] != "sample.name"){
               cat("Error: The column 1 of sample information must be sample.name")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               # cat("OK: There are no spaces in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "OK")
             }

             if(colnames(sample.info)[2] != "injection.order"){
               cat("Error: The column 2 of sample information must be injection.order")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               # cat("OK: There are no spaces in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "OK")
             }


             if(colnames(sample.info)[3] != "class"){
               cat("Error: The column 3 of sample information must be class")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               # cat("OK: There are no spaces in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "OK")
             }


             if(colnames(sample.info)[4] != "batch"){
               cat("Error: The column 4 of sample information must be batch")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               # cat("OK: There are no spaces in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "OK")
               if(any(table(sample.info[,"batch"]) <= 5)){
                 cat("There are some batch only have less than 5 samples. Please check it.")
                 sample.info.record <- c(sample.info.record, "Error")
               }
             }

             if(colnames(sample.info)[5] != "group"){
               cat("Error: The column 5 of sample information must be group")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               # cat("OK: There are no spaces in you sample.info.\n")
               sample.info.record <- c(sample.info.record, "OK")
             }

             class <-  unique(as.character(sample.info[,3]))
             if(all(c("QC", "Subject") %in% class)){
               sample.info.record <- c(sample.info.record, "OK")
             }else{
               cat("Error: The class of sample.information must be QC and Subject.\n")
               sample.info.record <- c(sample.info.record, "Error")
             }


             sample.idx <- match(sample.info$sample.name,
                                 unique(unlist(lapply(data, function(x) colnames(x)))))

             if(any(is.na(sample.idx))){
              cat("Error: The sample names in sample.inforamtion and data are not same.")
               sample.info.record <- c(sample.info.record, "Error")
             }else{
               sample.info.record <- c(sample.info.record, "OK")
             }


             cat("--------------------------------------------------------------\n")
             cat("Summary:\n")

             stat <- lapply(c(data.record, list(sample.info.record)),
                            function(x) {
                              return(c(ifelse(all(x!="Error"), "Valid", "Invalid"),
                                       sum(x=="OK"), sum(x=="Warning"), sum(x=="Error")))
                            })
             stat <- do.call(rbind, stat)
             colnames(stat) <- c("Check result","OK", "Warning", "Error")
             rownames(stat)[nrow(stat)] <- "sample.info"
             rownames(stat)[1:(nrow(stat)-1)] <- paste("batch", 1:(nrow(stat)-1), sep = "")

             print(stat, quote = FALSE)
             cat("\n")
             cat("\n")
             cat("data:\n")
             for(idx in 1:length(data.record)){
               if(all(data.record[idx] != "Error")){
                 cat(paste("Batch", idx), "is valid.\n")
               }else{
                 if(sum(data.record[idx] == "Warning") > 0){
                   cat("There", ifelse(sum(data.record[idx] == "Warning") > 1, "are", "is"),
                       sum(data.record[idx] == "Warning"),
                       ifelse(sum(data.record[idx] == "Warning") > 1, "Warnings", "Warning"),
                       "in your", paste("Batch", idx),
                       ".Please check it according to the information.\n")
                 }

                 if(sum(data.record[idx] == "Error") > 0){
                   cat("There", ifelse(sum(data.record[idx] == "Error") > 1, "are", "is"),
                       sum(data.record[idx] == "Error"),
                       ifelse(sum(data.record[idx] == "Error") > 1, "Errors", "Error"),
                       "in your", paste("Batch", idx),
                       ".Please check it according to the information.\n")
                 }

               }
             }



             cat("\n")
             cat("sample.info:\n")
             if(all(sample.info.record != "Error")){
               cat("sample.info is valid.\n")
             }else{
               if(sum(sample.info.record == "Warning") > 0){
                 cat("There", ifelse(sum(sample.info.record == "Warning") > 1, "are", "is"),
                     sum(sample.info.record == "Warning"),
                     ifelse(sum(sample.info.record == "Warning") > 1, "Warnings", "Warning"),
                     "in your sample.info. Please check it according to the information.\n")
               }

               if(sum(sample.info.record == "Error") > 0){
                 cat("There", ifelse(sum(sample.info.record == "Error") > 1, "are", "is"),
                     sum(sample.info.record == "Error"),
                     ifelse(sum(sample.info.record == "Error") > 1, "Errors", "Error"),
                     "in your sample.info. Please check it according to the information.\n")
               }

             }
             stat <- stat
           })





###version 0.9.73 backup
# setGeneric(name = "checkData",
#            def = function(data = "data.csv",
#                           sample.info = "sample.information.csv",
#                           path = "."){
#
#              data.record <- NULL
#              sample.info.record <- NULL
#
#              cat("Read data.\n")
#              cat("--------------------------------------------------------------\n")
#              data <- readr::read_csv(file.path(path, data), col_types = readr::cols(),
#                                      progress = TRUE)
#              data <- as.data.frame(data)
#              sample.info <- readr::read_csv(file.path(path, sample.info), col_types = readr::cols(),
#                                             progress = TRUE)
#
#              sample.info <- as.data.frame(sample.info)
#              #####
#
#              ##check data, component
#              data.col.name <- colnames(data)
#
#              if(all(data.col.name != "name")){
#                cat("Error: No name in data.\n")
#                data.record <- c(data.record, "Error")
#              }else{
#                # cat("OK: The first column of data is name.\n")
#                data.record <- c(data.record, "OK")
#              }
#
#              if(all(data.col.name != "mz")){
#                cat("Error: No mz in data.\n")
#                data.record <- c(data.record, "Error")
#              }else{
#                # cat("OK: The second column of data is mz.\n")
#                data.record <- c(data.record, "OK")
#              }
#
#              if(all(data.col.name != "rt")){
#                cat("Error: No rt in data.\n")
#                data.record <- c(data.record, "Error")
#              }else{
#                # cat("OK: The third column of data is not rt.\n")
#                data.record <- c(data.record, "OK")
#              }
#
#
#              ##check data, RT
#              # if(any(data.col.name == "rt")){
#              #   rt <- as.numeric(data$rt)
#              #   if(max(rt) < 60){
#              #     cat("Warning: Please confirm that the rt is measured in seconds,\n")
#              #     data.record <- c(data.record, "Warning")
#              #   }
#              #   cat("RT range is", range(rt), ".\n")
#              # }
#
#              cat("--------------------------------------------------------------\n")
#              ##check sample.info
#              # if(ncol(sample.info) > 5){
#              #   cat("Error: The smple.info has more than two columns. Please check it.\n")
#              #   sample.info.record <- c(sample.info.record, "Error")
#              # }
#
#              if(sum(is.na(sample.info)) > 0){
#                cat("Error: There are", sum(is.na(sample.info)), "NAs in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                # cat("OK: There are no NAs in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "OK")
#              }
#
#              if(ifelse(is.na(sum(sample.info == "") > 0), FALSE, sum(sample.info == "") > 0)){
#                cat("Error: There are", sum(sample.info == ""), "spaces in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                # cat("OK: There are no spaces in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "OK")
#              }
#
#
#              if(colnames(sample.info)[1] != "sample.name"){
#                cat("Error: The column 1 of sample information must be sample.name")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                # cat("OK: There are no spaces in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "OK")
#              }
#
#              if(colnames(sample.info)[2] != "injection.order"){
#                cat("Error: The column 2 of sample information must be injection.order")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                # cat("OK: There are no spaces in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "OK")
#              }
#
#
#              if(colnames(sample.info)[3] != "class"){
#                cat("Error: The column 3 of sample information must be class")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                # cat("OK: There are no spaces in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "OK")
#              }
#
#
#              if(colnames(sample.info)[4] != "batch"){
#                cat("Error: The column 4 of sample information must be batch")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                # cat("OK: There are no spaces in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "OK")
#                if(any(table(sample.info[,"batch"]) <= 5)){
#                  cat("There are some batch only have less than 5 samples. Please check it.")
#                  sample.info.record <- c(sample.info.record, "Error")
#                }
#              }
#
#              if(colnames(sample.info)[5] != "group"){
#                cat("Error: The column 5 of sample information must be group")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                # cat("OK: There are no spaces in you sample.info.\n")
#                sample.info.record <- c(sample.info.record, "OK")
#              }
#
#
#              # group <-  table(as.character(sample.info[,5]))
#              # group1 <- matrix(group, nrow = 1)
#              # colnames(group1) <- names(group)
#              # rownames(group1) <- "Number"
#              #
#              # cat("Group information:\n")
#              # print(group1)
#
#
#              class <-  unique(as.character(sample.info[,3]))
#              if(all(c("QC", "Subject") %in% class)){
#                sample.info.record <- c(sample.info.record, "OK")
#              }else{
#                cat("Error: The class of sample.information must be QC and Subject.\n")
#                sample.info.record <- c(sample.info.record, "Error")
#              }
#
#              # if(any(group1[1,] < 2)){
#              #   cat("Error: ",group1[1,which(group1[,1] < 3)],
#              #       ifelse(length(which(group1[,1]) < 3) > 1, "have", "has"),
#              #       "less than 3 samples, please check it.\n"
#              #   )
#              #   sample.info.record <- c(sample.info.record, "Error")
#              # }else{
#              #   sample.info.record <- c(sample.info.record, "OK")
#              # }
#
#
#              sample.idx <- match(sample.info$sample.name, colnames(data))
#
#              if(any(is.na(sample.idx))){
#                cat("Error: The sample names in sample.inforamtion and data are not same.")
#                sample.info.record <- c(sample.info.record, "Error")
#              }else{
#                sample.info.record <- c(sample.info.record, "OK")
#              }
#
#              # data.col.name <- setdiff(data.col.name, c("name", "mz", "rt"))
#              # sample.name <- as.character(sample.info[,1])
#              #
#              # if(any(data.col.name != sample.name)){
#              #   cat("Error: The sample names in data and sample.info are not same.\n")
#              #   sample.info.record <- c(sample.info.record, "Error")
#              # }else{
#              #   sample.info.record <- c(sample.info.record, "OK")
#              # }
#
#              # cat("--------------------------------------------------------------\n")
#              ##ms2.file
#              # ms2.file <- grep(ifelse(ms2.type == "mgf","mgf", "mzXML"), dir(path), value = TRUE)
#              # if(length(ms2.file) == 0){
#              #   cat("Error: There are no ms2 data.\n")
#              #   ms2.file.record <- c(ms2.file.record, "Error")
#              # }else{
#              #   cat("There are", length(ms2.file), "ms2 data.\n")
#              #   ms2.file.record <- c(ms2.file.record, "OK")
#              #   cat("Total size of ms2 data is", sum(file.size(ms2.file)/(1024*1024)), "M.\n")
#              # }
#
#              cat("--------------------------------------------------------------\n")
#              cat("Summary:\n")
#
#              stat <- lapply(list(data.record, sample.info.record),
#                             function(x) {
#                               return(c(ifelse(all(x!="Error"), "Valid", "Invalid"),
#                                        sum(x=="OK"), sum(x=="Warning"), sum(x=="Error")))
#                             })
#              stat <- do.call(rbind, stat)
#              colnames(stat) <- c("Check result","OK", "Warning", "Error")
#              rownames(stat) <- c("data", "sample.info")
#
#              print(stat, quote = FALSE)
#              cat("\n")
#
#              cat("\n")
#              cat("data:\n")
#              if(all(data.record != "Error")){
#                cat("data is valid.\n")
#              }else{
#                if(sum(data.record == "Warning") > 0){
#                  cat("There", ifelse(sum(data.record == "Warning") > 1, "are", "is"),
#                      sum(data.record == "Warning"),
#                      ifelse(sum(data.record == "Warning") > 1, "Warnings", "Warning"),
#                      "in your data. Please check it according to the information.\n")
#                }
#
#                if(sum(data.record == "Error") > 0){
#                  cat("There", ifelse(sum(data.record == "Error") > 1, "are", "is"),
#                      sum(data.record == "Error"),
#                      ifelse(sum(data.record == "Error") > 1, "Errors", "Error"),
#                      "in your data. Please check it according to the information.\n")
#                }
#
#              }
#
#              cat("\n")
#              cat("sample.info:\n")
#              if(all(sample.info.record != "Error")){
#                cat("sample.info is valid.\n")
#              }else{
#                if(sum(sample.info.record == "Warning") > 0){
#                  cat("There", ifelse(sum(sample.info.record == "Warning") > 1, "are", "is"),
#                      sum(sample.info.record == "Warning"),
#                      ifelse(sum(sample.info.record == "Warning") > 1, "Warnings", "Warning"),
#                      "in your sample.info. Please check it according to the information.\n")
#                }
#
#                if(sum(sample.info.record == "Error") > 0){
#                  cat("There", ifelse(sum(sample.info.record == "Error") > 1, "are", "is"),
#                      sum(sample.info.record == "Error"),
#                      ifelse(sum(sample.info.record == "Error") > 1, "Errors", "Error"),
#                      "in your sample.info. Please check it according to the information.\n")
#                }
#
#              }
#              stat <- stat
#            })

# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyBS)
library(cttools)
library(extrantsr)
library(fslr)
library(matlabr)
if (!have_matlab()) {
  options(matlab.path = '/Applications/MATLAB_R2014b.app/bin')
}

# filename =  "~/Desktop/6003-157/_AXIAL_HEAD_STD_20140115181000_2_Eq_1.nii.gz"
verbose = TRUE

# Simple Function that will get the messages out
add_expr = function(expr){
  withCallingHandlers(
    expr,
    message = function(m) {
      shinyjs::text(id = "nText", 
                    text = paste0(m$message, "<br/>"), add = TRUE)
    },
    warning = function(w) {
      shinyjs::text(id = "nText", 
                    text = paste0(w$message, "<br/>"), add = TRUE)
    },
    error = function(e) {
      shinyjs::text(id = "nText", 
                    text = paste0(e$message, "<br/>"), add = TRUE)
    })
}

shinyServer(function(input, output, session) {
  values <- reactiveValues(ss = NULL,
                           preprocess = NULL, 
                           df = NULL,
                           nim = NULL,
                           res = NULL)
  
  
  read_img = reactive({
    inFile <- input$img_fname
    
    if (is.null(inFile)) {
      return(NULL)
    } else {
      fname = inFile$name
      
      dpath = inFile$datapath
      #######################
      # Renaming Images so they have the filename and .nii.gz
      #######################        
      filename = file.path(dirname(dpath), fname)
      if (file.exists(dpath)) {
        if (file.exists(filename)) {
          file.remove(filename)
        }
        file.rename(dpath, filename)
      }
      add_expr({
        if (verbose) {
          message("<h3># Reading in Data</h3>")
        }       
        img = readnii(filename)
      })
      return(img)
    }
  })
  #   
  #   proc_img = reactive({
  #     
  #     img = read_img()
  #     if (!is.null(img)) {
  #       preprocess = ich_preprocess(img = img, 
  #                                   mask = NULL,
  #                                   verbose = TRUE)
  #       return(preprocess)
  #     } else {
  #       return(NULL)
  #     }
  #     
  #   })
  
  output$origPlot <- renderPlot({
    shinyjs::disable("download")
    
    img = read_img()
    if (!is.null(img)) {
      add_expr({
        if (verbose) {
          message("<h3># Plotting Original in Data</h3>")
        }     
      })        
      ortho2(img, window = c(0, 100), text = "Original\nImage")
      
      updateButton(session, "ss", 
                   disabled = FALSE, 
                   icon = icon("chevron-right"))
    } else {
      return()
    }
  })
  
  #########################################
  # Skull Stripping
  #########################################  
  observe({
    
    if (input$ss > 0) {
      
      updateButton(session, "ss", disabled = TRUE, 
                   icon = icon("spinner fa-spin"))      
      
      img = read_img()
      add_expr({
        if (verbose) {
          message("<h3># Skull-Stripping Images</h3>")
        } 
        ss = CT_Skull_Strip_robust(img, retimg = TRUE)
        # ss = CT_Skull_Strip(img, retimg = TRUE)
      })
      mask = ss > 0
      img = check_nifti(img)
      ss = mask_img(img, mask)
      ss = window_img(ss, window = c(0, 100), replace = "zero")
      values$ss = ss
      icon = icon("check")
      updateButton(session, "ss", 
                   disabled = TRUE, 
                   icon = icon("check"))          
      updateButton(session, "reg", disabled = FALSE, 
                   icon = icon("chevron-right"))      
    }
  })
  
  #########################################
  # Skull Stripping Plot
  #########################################    
  output$ssPlot <- renderPlot({
    # print(input$ss)
    ss = values$ss
    if (!is.null(ss)) {
      add_expr({
        if (verbose) {
          message("<h3># Plotting Skull-Stripped Data<h3>")
        }    
      })        
      ortho2(ss, window = c(0, 100), xyz = xyz(ss > 0), 
             text = "Skull-Stripped\nImage")
      
    } else {
      return()
    }
  })  
  
  
  #########################################
  # Registration
  #########################################    
  observe({
    # print(input$reg)
    if (input$reg > 0) {    
      updateButton(session, "reg", disabled = TRUE, 
                   icon = icon("spinner fa-spin"))      
      
      ss = values$ss
      mask = ss > 0
      outprefix = tempfile()
      typeofTransform = "Rigid"
      template.file = system.file(
        "scct_unsmooth_SS_0.01.nii.gz", 
        package = "cttools")
      template = readnii(template.file)
      interpolator = "Linear"
      add_expr({
        if (verbose) {
          message("<h3># Registering Images<h3>")
        }            
        preprocess = registration(
          filename = ss, 
          skull_strip = FALSE,
          correct = FALSE, 
          outfile = NULL, 
          retimg = TRUE, 
          typeofTransform = typeofTransform,
          template.file = template.file, 
          interpolator = interpolator,
          remove.warp = FALSE,
          outprefix = outprefix,
          verbose = verbose)
        
        omask = ants_apply_transforms(
          fixed = template.file,
          moving = mask, 
          typeofTransform = typeofTransform,
          interpolator = interpolator, 
          transformlist = preprocess$fwdtransforms)
      })
      preprocess$mask = mask
      preprocess$ss_image = ss
      preprocess$transformed_image = preprocess$outfile
      preprocess$transformed_mask = omask
      preprocess$outfile = NULL
      preprocess$fixed = template
      
      values$preprocess = preprocess
      updateButton(session, "reg", disabled = TRUE, 
                   icon = icon("check"))
      updateButton(session, "make_pred", disabled = FALSE, 
                   icon = icon("chevron-right"))      
      
    }
  })
  
  
  #########################################
  # Registration Plot
  #########################################      
  output$regPlot <- renderPlot({
    # print(input$ss)
    preprocess = values$preprocess
    if (!is.null(preprocess)) {
      add_expr({
        if (verbose) {
          message("<h3># Plotting Registered Data<h3>")
        }      
        double_ortho(preprocess$transformed_image, 
                     preprocess$fixed, 
                     window = c(0, 100), 
                     text = "Registered\nImage\nand\nTemplate")
      })
    } else {
      return()
    }
  })
  
  
  #########################################
  # Making Predictors
  #########################################      
  observe({
    print(input$make_pred)
    if (input$make_pred > 0) {    
      preprocess = values$preprocess
      updateButton(session, "make_pred", disabled = TRUE, 
                   icon = icon("spinner fa-spin"))
      
      timg = preprocess$transformed_image
      tmask = preprocess$transformed_mask > 0.5
      add_expr({
        if (verbose) {
          message("<h3># Making Predictors<h3>")
        }         
        img.pred = make_predictors(timg, 
                                   mask = tmask,
                                   roi = NULL, 
                                   save_imgs = FALSE,
                                   verbose = verbose)
      })
      df = img.pred$df
      nim = img.pred$nim
      rm(list = "img.pred")
      gc()
      df$multiplier = ich_candidate_voxels(df)

      
      values$df = df
      values$nim = nim
      updateButton(session, "make_pred", 
                   disabled = TRUE, icon = icon("check"))
      updateButton(session, "predict", disabled = FALSE, 
                   icon = icon("chevron-right"))      
    }
  })
  
  #########################################
  # Making Predictors Plot
  #########################################          
  output$candPlot <- renderPlot({
    
    df = values$df
    nim = values$nim    
    preprocess = values$preprocess
    
    if (!is.null(preprocess) & !is.null(df) & !is.null(nim)) {
      cand = remake_img(df$multiplier, nim)
      xyz = xyz(cand)
      timg = preprocess$transformed_image
      add_expr({
        if (verbose) {
          message("<h3># Plotting Candidate Data<h3>")
        }     
      })
      ortho2(timg, cand, 
             window = c(0, 100), 
             text = "Registered\nImage\nand\nCandidate Mask")
    } else {
      return()
    }
  })  
  
  #########################################
  # Predict Image
  #########################################      
  observe({
    print(input$predict)
    if (input$predict > 0) {  
      img = read_img()
      df = values$df
      nim = values$nim
      preprocess = values$preprocess
      updateButton(session, "predict", disabled = TRUE, 
                   icon = icon("spinner fa-spin"))
      
      model = "rf"
      add_expr({
        
        if (verbose) {
          message("<h3># Making Prediction Images<h3>")
        }    
        L = ich_predict(df = df, 
                        native_img = img, 
                        nim = nim,
                        model = model,
                        verbose = verbose,
                        transformlist = preprocess$invtransforms,
                        interpolator = preprocess$interpolator)
      })
      L$preprocess = preprocess
      values$res = L
      updateButton(session, "predict", disabled = TRUE, 
                   icon = icon("check"))
      shinyjs::enable("download")
    }
  })
  
  #########################################
  # Predict Image Plot
  #########################################          
  output$predPlot <- renderPlot({
    res = values$res
    scc = res$native_prediction$smoothed_prediction_image
    # cc = res$native_prediction$prediction_image
    ss = values$ss
    
    if (!is.null(ss) & !is.null(scc)) {
      xyz = xyz(scc > 0.5)
      add_expr({
        if (verbose) {
          message("<h3># Plotting Prediction Mask<h3>")
        }     
      })         
      ortho2(ss, scc, 
             window = c(0, 100),
             text = "Skull-Stripped\nImage\nand\nPrediction Mask")
    } else {
      return()
    }
  })    
  
  output$download <- downloadHandler(
    filename = 'output.rda',
    content = function(file) {
      isolate({
        vv = reactiveValuesToList(values)
      })
      res = vv$res
      
      save(res, file = file)
    }
  )
  
  output$download_pred <- downloadHandler(
    filename = 'smoothed_prediction_image.nii.gz',
    content = function(file) {
      isolate({
        vv = reactiveValuesToList(values)
      })
      res = vv$res
      
      save(res, file = file)
    }
  )  
  
  output$download_pred <- downloadHandler(
    filename = 'smoothed_prediction_image.nii.gz',
    content = function(file) {
      isolate({
        vv = reactiveValuesToList(values)
      })
      outimg = vv$res$native_prediction$smoothed_prediction_image
      writenii(outimg, 
               filename = file, 
               gzipped = TRUE)
      file.rename(paste0(file, ".nii.gz"), file)
    }
  )

  
})



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

remake_img = function(vec, img, mask = NULL) {
  if (is.null(mask)) {
    mask = array(1, dim = dim(img))
  }
  img2 = niftiarr(img, 0)
  img2[mask == 1] = vec
  img2 = datatyper(img2)
  img2 = cal_img(img2)
  return(img2)
}

shinyServer(function(input, output, session) {
  values <- reactiveValues(ss = NULL,
                           preprocess = NULL, 
                           df = NULL,
                           nim = NULL,
                           res = NULL)
  
  
  read_img = reactive({
    inFile <- input$img_fname
    
    # print(inFile)
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
      img = readnii(filename)
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
    print(input$ss)
    if (input$ss > 0) {
      
      updateButton(session, "ss", disabled = TRUE, 
                   icon = icon("spinner fa-spin"))      
      
      img = read_img()
      
      ss = CT_Skull_Strip_robust(img, retimg = TRUE)
      mask = ss > 0
      img = check_nifti(img)
      ss = mask_img(img, mask)
      ss = window_img(ss, window = c(0, 100), replace = "zero")
      values$ss = ss
      icon = icon("check")
      updateButton(session, "ss", disabled = disabled, icon = icon("check"))          
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
    print(input$reg)
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
      verbose = TRUE
      
      res = registration(
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
      
      omask = ants_apply_transforms(fixed = template.file,
                                    moving = mask, 
                                    typeofTransform = typeofTransform,
                                    interpolator = interpolator, 
                                    transformlist = res$fwdtransforms)
      res$mask = mask
      res$ss_image = ss
      res$transformed_image = res$outfile
      res$transformed_mask = omask
      res$outfile = NULL
      res$fixed = template
      
      values$preprocess = res
      updateButton(session, "reg", disabled = TRUE, icon = icon("check"))
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
      double_ortho(preprocess$transformed_image, 
             preprocess$fixed, 
             window = c(0, 100), 
             text = "Registered\nImage\nand\nTemplate")
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
      verbose = TRUE
      img.pred = make_predictors(timg, 
                                 mask = tmask,
                                 roi = NULL, 
                                 save_imgs = FALSE,
                                 verbose = verbose)
      df = img.pred$df
      nim = img.pred$nim
      rm(list = "img.pred")
      gc()
      df$multiplier = ich_candidate_voxels(df)
      values$df = df
      values$nim = nim
      updateButton(session, "make_pred", disabled = TRUE, icon = icon("check"))
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
      verbose = TRUE
      if (verbose) {
        message("# Making Prediction Images")
      }    
      L = ich_predict(df = df, 
                      nim = nim, 
                      model = model,
                      verbose = verbose,
                      transformlist = preprocess$invtransforms,
                      interpolator = preprocess$invtransforms)
      L$preprocess = preprocess
      values$res = L
      updateButton(session, "predict", disabled = TRUE, icon = icon("check"))
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
      res = values$res
      save(res, file = file)
    }
  )
  
  #   proc_img = reactive({
  #     
  #     img = read_img()
  #     if (!is.null(img)) {
  #       preprocess = ich_preprocess(img = img, 
  #                                mask = NULL,
  #                                verbose = TRUE)
  #       return(preprocess)
  #     } else {
  #       return(NULL)
  #     }
  #     
  #   })    
  
  #   output$tout <- renderPlot({
  #     proc_img()
  #   })
  
  
  
})


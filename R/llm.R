# LLM Provider System - Single file for all LLM functionality

# Internal: retry a function call with exponential backoff on HTTP 429 errors.
# fn          - a zero-argument function that makes the API call
# max_retries - maximum number of retry attempts (not counting the first try)
# initial_wait - seconds to wait before the first retry (doubles each attempt)
.retry_with_backoff <- function(fn, max_retries = 3, initial_wait = 5) {
  attempt <- 0L
  repeat {
    result <- tryCatch(
      list(value = fn(), error = NULL),
      error = function(e) list(value = NULL, error = e)
    )
    if (is.null(result$error)) {
      return(result$value)
    }
    err_msg <- conditionMessage(result$error)
    is_rate_limit <- grepl("429|Too Many Requests|rate.?limit", err_msg, ignore.case = TRUE)
    attempt <- attempt + 1L
    if (!is_rate_limit || attempt > max_retries) {
      stop(result$error)
    }
    wait_secs <- initial_wait * (2L ^ (attempt - 1L))
    message(sprintf(
      "[OntoAnno] Rate limit (429). Waiting %d s before retry %d/%d ...",
      wait_secs, attempt, max_retries
    ))
    Sys.sleep(wait_secs)
  }
}

# Internal: summarize params object for diagnostics.
.summarize_llm_params <- function(params) {
  if (is.null(params)) {
    return("none")
  }
  param_names <- tryCatch(names(params), error = function(e) NULL)
  if (is.null(param_names) || length(param_names) == 0) {
    return("provided (names unavailable)")
  }
  paste(param_names, collapse = ", ")
}

#' Call LLM API
#'
#' Routes prompts to OpenAI, Anthropic, Google Gemini, or local Ollama models via ellmer package.
#'
#' @param prompt Character. The prompt to send.
#' @param provider Character. One of "openai", "anthropic", "gemini", "ollama", "vllm",
#'   "claude_cli" (Claude Code CLI in headless `-p` mode) or "codex_cli" (Codex CLI `exec` mode).
#'   The two CLI providers need no API key or HTTP endpoint; they call the locally installed
#'   official CLI, which uses its own login.
#' @param model Character. Model name (e.g., "gpt-5.2", "claude-opus-4-6", "gemini-2.0-flash").
#'   For Ollama, use model names like "llama2", "mistral", "neural-chat", etc.
#' @param params An ellmer params object created by ellmer::params().
#' @param api_key Character. API key (uses env var if NULL).
#'   OpenAI: OPENAI_API_KEY
#'   Anthropic: ANTHROPIC_API_KEY
#'   Gemini: GOOGLE_API_KEY or GEMINI_API_KEY
#'   Ollama/vLLM: Not required (local)
#' @param system_prompt Character. System prompt.
#' @param api_url Character. API URL for local LLMs (Ollama or vLLM). For the CLI providers,
#'   the path to the CLI binary (default: "claude" or "codex" found on PATH).
#'   Ollama default: "http://localhost:11434"
#'   vLLM default: "http://localhost:8000"
#'
#' @details
#' **Ollama Setup:**
#' To use Ollama locally, first install it from https://ollama.ai, then run:
#' ```
#' ollama pull llama2       # or another model
#' ollama serve             # start the server (runs on localhost:11434)
#' ```
#' Popular models: llama2, mistral, neural-chat, dolphin-mixtral, etc.
#'
#' **vLLM Setup:**
#' Similar to Ollama but optimized for performance. Requires Python installation.
#'
#' @return Character. The LLM response text.
#' @importFrom ellmer params chat_openai chat_anthropic chat_google_gemini chat_ollama chat_vllm
#' @export
call_llm <- function(prompt, provider = "openai", model = NULL,
                     params = NULL, api_key = NULL, system_prompt = NULL,
                     api_url = NULL) {
  # Ollama models often include a tag suffix (e.g. "llama3.2:latest").
  # The ellmer client expects the bare model name, so strip any tag.
  if (!is.null(model) && provider == "ollama") {
    model <- sub(":.*$", "", model)
  }
  if (!provider %in% c("claude_cli", "codex_cli") && !requireNamespace("ellmer", quietly = TRUE)) {
    stop("Package 'ellmer' required but not installed")
  }
  # warn if provider uses API key and none provided
  if (provider %in% c("openai", "anthropic", "gemini")) {
    envs <- switch(provider,
                   openai = "OPENAI_API_KEY",
                   anthropic = "ANTHROPIC_API_KEY",
                   gemini = c("GOOGLE_API_KEY", "GEMINI_API_KEY"))
    if (is.null(api_key)) {
      has <- vapply(envs, function(e) nzchar(Sys.getenv(e, "")), logical(1))
      if (!any(has)) {
        message("[OntoAnno] warning: no API key found for provider '", provider, "'. ",
                "Set one via environment variable (", paste(envs, collapse = "/"), ") ",
                "or pass via llm_config$api_key.")
      }
    }
  }

  if (provider == "openai") {
    chat <- ellmer::chat_openai(
      system_prompt = system_prompt,
      model = model,
      params = params,
      api_key = api_key,
      echo = "none"
    )
  } else if (provider == "anthropic") {
    chat <- ellmer::chat_anthropic(
      system_prompt = system_prompt,
      model = model,
      params = params,
      api_key = api_key,
      echo = "none"
    )
  } else if (provider == "gemini") {
    chat <- ellmer::chat_google_gemini(
      system_prompt = system_prompt,
      model = model,
      params = params,
      api_key = api_key,
      echo = "none"
    )
  } else if (provider == "ollama") {
    # Set default Ollama API URL if not provided
    if (is.null(api_url)) {
      api_url <- "http://localhost:11434"
    }
    chat <- ellmer::chat_ollama(
      system_prompt = system_prompt,
      model = model,
      params = params,
      base_url = api_url,
      echo = "none"
    )
  } else if (provider == "vllm") {
    # Set default vLLM API URL if not provided
    if (is.null(api_url)) {
      api_url <- "http://localhost:8000"
    }
    chat <- ellmer::chat_vllm(
      system_prompt = system_prompt,
      model = model,
      params = params,
      base_url = api_url,
      echo = "none"
    )
  } else if (provider %in% c("claude_cli", "codex_cli")) {
    # Official headless modes of the Claude Code and Codex CLIs: no HTTP API, no key;
    # the CLI's own login is used. `params` are ignored (the CLIs expose no sampling
    # controls); `system_prompt` is passed via --append-system-prompt (claude) or
    # prepended to the prompt (codex). The chat object mimics ellmer's `$chat()`.
    bin <- if (!is.null(api_url)) api_url else if (provider == "claude_cli") "claude" else "codex"
    if (!nzchar(Sys.which(bin)) && !file.exists(bin)) {
      stop("CLI binary not found: ", bin, ". Pass its path via llm_config$api_url.")
    }
    chat <- list(chat = function(p) {
      if (provider == "claude_cli") {
        args <- c("-p", "--output-format", "text", "--no-session-persistence")
        if (!is.null(model)) args <- c(args, "--model", model)
        if (!is.null(system_prompt)) args <- c(args, "--append-system-prompt", shQuote(system_prompt))
        out <- withCallingHandlers(
          system2(bin, args, input = p, stdout = TRUE, stderr = TRUE),
          warning = function(w) invokeRestart("muffleWarning"))
      } else {
        outfile <- tempfile(fileext = ".txt")
        on.exit(unlink(outfile), add = TRUE)
        args <- c("exec", "--skip-git-repo-check", "-s", "read-only", "-o", outfile)
        if (!is.null(model)) args <- c(args, "-m", model)
        full <- if (is.null(system_prompt)) p else paste(system_prompt, p, sep = "\n\n")
        withCallingHandlers(
          system2(bin, args, input = full, stdout = FALSE, stderr = FALSE),
          warning = function(w) invokeRestart("muffleWarning"))
        out <- if (file.exists(outfile)) readLines(outfile, warn = FALSE) else character()
      }
      status <- attr(out, "status")
      if (!is.null(status) && status != 0) stop(provider, " exited with status ", status, ": ", paste(out, collapse = "\n"))
      if (!length(out)) stop(provider, " returned no output")
      paste(out, collapse = "\n")
    })
  } else {
    stop("Unknown provider: ", provider, ". Use 'openai', 'anthropic', 'gemini', 'ollama', 'vllm', 'claude_cli' or 'codex_cli'")
  }

  tryCatch({
    .retry_with_backoff(function() chat$chat(prompt))
  }, error = function(e) {
    err_msg <- conditionMessage(e)
    params_used <- .summarize_llm_params(params)

    msg <- paste0(
      "LLM request failed (provider=", provider,
      ", model=", model %||% "<default>",
      ", params=", params_used, ").\n",
      err_msg
    )

    if (identical(provider, "openai") && grepl("HTTP 400", err_msg, fixed = TRUE)) {
      msg <- paste0(
        msg,
        "\nHint: This often indicates one or more params are unsupported for the selected model.",
        " Check model docs and retry with fewer params."
      )
    }

    stop(msg, call. = FALSE)
  })
}

#' Prepare LLM Configuration
#'
#' Helper function that merges legacy `model` parameter with new `llm_config` structure.
#' Provides backward compatibility while supporting new multi-provider configuration.
#'
#' @param model Character. Model name (default: "gpt-5.2"). Used if llm_config not provided.
#' @param llm_config Optional list. Configuration with: provider, model, params,
#'   temperature, max_tokens, api_key, api_url, system_prompt.
#'   If `params` is provided, it should be an `ellmer::params(...)` object and is
#'   passed directly to the chat client.
#'
#' @return A list with merged LLM configuration.
#'
#' @keywords internal
#' @export
prepare_config <- function(model = "gpt-5.2", llm_config = NULL) {
  config <- list(
    provider = "openai",
    model = model,
    params = NULL,
    temperature = NULL,
    max_tokens = NULL,
    api_key = NULL,
    api_url = NULL,
    system_prompt = NULL
  )

  if (!is.null(llm_config)) {
    for (key in names(llm_config)) {
      config[[key]] <- llm_config[[key]]
    }
  }

  return(config)
}

# Helper: NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Check Ollama Service Availability
#'
#' Checks if Ollama is running and accessible at the specified URL.
#'
#' @param api_url Character. Ollama API URL (default: "http://localhost:11434")
#'
#' @return Logical. TRUE if Ollama is accessible, FALSE otherwise.
#'
#' @details
#' This function attempts to ping the Ollama API to verify it is running.
#' Returns FALSE silently if the service is not available.
#'
#' @examples
#' \dontrun{
#' check_ollama_available()
#' }
#'
#' @export
check_ollama_available <- function(api_url = "http://localhost:11434") {
  tryCatch({
    response <- httr::GET(paste0(api_url, "/api/tags"), httr::timeout(2))
    httr::status_code(response) == 200
  }, error = function(e) {
    FALSE
  })
}

#' List Available Ollama Models
#'
#' Lists all models currently available in Ollama.
#'
#' @param api_url Character. Ollama API URL (default: "http://localhost:11434")
#'
#' @return Character vector of model names, or NULL if Ollama is not available.
#'
#' @details
#' Requires Ollama to be running. If Ollama is not accessible, returns NULL with a warning.
#'
#' @examples
#' \dontrun{
#' list_ollama_models()
#' # Returns something like:
#' # [1] "llama2" "mistral" "neural-chat"
#' }
#'
#' @export
list_ollama_models <- function(api_url = "http://localhost:11434") {
  if (!check_ollama_available(api_url)) {
    warning("Ollama is not running at ", api_url, ". ",
            "Start Ollama with: ollama serve")
    return(NULL)
  }
  
  tryCatch({
    response <- httr::GET(paste0(api_url, "/api/tags"))
    content <- httr::content(response, as = "parsed")
    if (!is.null(content$models)) {
      return(sapply(content$models, function(x) x$name))
    }
    return(NULL)
  }, error = function(e) {
    warning("Error fetching Ollama models: ", e$message)
    return(NULL)
  })
}

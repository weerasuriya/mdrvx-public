  }

  if (mode == 0) {
    # Local mode - return detailed output
    return(output)
  } else if (mode == 1) {
    # Cluster mode - return only fitting values
    return(fo)
  }

}

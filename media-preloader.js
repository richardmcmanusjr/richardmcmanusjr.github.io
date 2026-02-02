// Comprehensive Media Preloader - Preloads images and videos for both main and running sites
(function() {
  // Configuration for media preloading
  const config = {
    // Main site media
    mainSite: {
      images: [
        'hero-image.png',
        'AirSpacerTestStructure2.png',
        'radiation.png',
        'McManus_Richard_CV.pdf',
        'GOATCoach.png',
        'Charles.png',
        'coach.png',
        'berkeleyHalf.png'
      ],
      videos: [
        'running.mov',
        'hiking.mov',
        'weightlifting.mov',
        'efoil.mov'
      ]
    },
    // Running site media
    runningSite: {
      images: [
        'GOATCoach.png',
        'Charles.png',
        'coach.png',
        'berkeleyHalf.png'
      ],
      videos: [
        'video1.mov',
        'video2.mov',
        'video3.mov'
      ]
    }
  };

  // Preload settings
  const preloadSettings = {
    imageBatchSize: 5,
    videoBatchSize: 2,
    batchDelay: 100, // ms between batches
    videoPreloadTimeout: 5000, // ms before forcing video to be considered loaded
    showProgress: true // Set to true for console logging
  };

  // Utility function to detect which site we're on
  function detectCurrentSite() {
    const currentPath = window.location.pathname.toLowerCase();
    if (currentPath.includes('/running.html') || currentPath.includes('/running/')) {
      return 'running';
    }
    return 'main';
  }

  // Get media list based on current site
  function getMediaList() {
    const site = detectCurrentSite();
    return config[site === 'running' ? 'runningSite' : 'mainSite'];
  }

  // Preload images
  function preloadImages(imageList, callback) {
    if (!imageList || imageList.length === 0) {
      if (callback) callback();
      return;
    }

    let loadedCount = 0;
    let currentBatch = 0;

    function loadBatch() {
      const start = currentBatch * preloadSettings.imageBatchSize;
      const end = Math.min(start + preloadSettings.imageBatchSize, imageList.length);

      for (let i = start; i < end; i++) {
        const img = new Image();
        img.onload = img.onerror = () => {
          loadedCount++;
          if (preloadSettings.showProgress) {
            console.log(`Image loaded: ${imageList[i]} (${loadedCount}/${imageList.length})`);
          }
          if (loadedCount === imageList.length && callback) {
            callback();
          }
        };
        img.src = imageList[i];
      }

      currentBatch++;

      if (end < imageList.length) {
        setTimeout(loadBatch, preloadSettings.batchDelay);
      }
    }

    loadBatch();
  }

  // Preload videos
  function preloadVideos(videoList, callback) {
    if (!videoList || videoList.length === 0) {
      if (callback) callback();
      return;
    }

    let loadedCount = 0;
    let currentBatch = 0;

    function loadBatch() {
      const start = currentBatch * preloadSettings.videoBatchSize;
      const end = Math.min(start + preloadSettings.videoBatchSize, videoList.length);

      for (let i = start; i < end; i++) {
        const videoPath = videoList[i];
        const video = document.createElement('video');
        video.preload = 'auto';
        video.style.display = 'none';

        let videoLoaded = false;

        function markLoaded() {
          if (!videoLoaded) {
            videoLoaded = true;
            loadedCount++;
            if (preloadSettings.showProgress) {
              console.log(`Video preloaded: ${videoPath} (${loadedCount}/${videoList.length})`);
            }
            if (loadedCount === videoList.length && callback) {
              callback();
            }
          }
        }

        // Set up event listeners
        video.addEventListener('canplay', markLoaded, { once: true });
        video.addEventListener('error', markLoaded, { once: true });

        // Timeout fallback
        setTimeout(markLoaded, preloadSettings.videoPreloadTimeout);

        // Add source and start loading
        const source = document.createElement('source');
        source.src = videoPath;
        source.type = 'video/mp4';
        video.appendChild(source);

        // Start loading
        video.load();
      }

      currentBatch++;

      if (end < videoList.length) {
        setTimeout(loadBatch, preloadSettings.batchDelay);
      }
    }

    loadBatch();
  }

  // Main preload function
  function preloadAllMedia() {
    const mediaList = getMediaList();
    const site = detectCurrentSite();

    if (preloadSettings.showProgress) {
      console.log(`Starting media preload for ${site} site...`);
    }

    // Preload images first, then videos
    preloadImages(mediaList.images, () => {
      if (preloadSettings.showProgress) {
        console.log('All images preloaded');
      }
      preloadVideos(mediaList.videos, () => {
        if (preloadSettings.showProgress) {
          console.log('All media preloaded successfully');
        }
        // Dispatch custom event to signal completion
        window.dispatchEvent(new CustomEvent('mediaPreloadComplete'));
      });
    });
  }

  // Start preloading when DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', preloadAllMedia);
  } else {
    preloadAllMedia();
  }

  // Also trigger on window load as fallback
  window.addEventListener('load', preloadAllMedia);

  // Expose utility functions to window for debugging if needed
  window.mediaPreloader = {
    preloadImages,
    preloadVideos,
    preloadAllMedia,
    getMediaList,
    getConfig: () => config,
    getSettings: () => preloadSettings
  };
})();

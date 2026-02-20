// Image Preloader - Preloads images for smoother carousel and page transitions
(function() {
  // All images to preload
  const imagesToPreload = [
    // Team images
    'content/team/salahuddin.jpg',
    'content/team/KoushikDas.jpg',
    'content/team/ChiragGarg.jpg',
    'content/team/Chia-ChunLee.png',
    'content/team/RamitDutta.jpg',
    'content/team/RichardMcManus.png',
    'content/team/KimPham.jpg',
    'content/team/AniketSadishiva.jpg',
    'content/team/Wei-YangWeng.jpg',
    'content/team/YonghyeokSon.jpg',
    'content/team/JooseongYun.jpg',
    'content/team/You.png',
    'content/team/JunghyeonHwang.jpg',
    'content/team/MyeongseopSong.jpg',
    'content/team/UrmitaSikder.png',
    'content/team/RakshitJain.png',
    'content/team/SteveVolkman.png',
    // Research icons
    'content/icon-materials.png',
    'content/icon-simulation.png',
    'content/icon-engineering.png',
    'content/icon-engineering2.png',
    'content/icon-computing.png',
    // Activity images
    'content/Food1.png',
    'content/Food2.png',
    'content/Food3.png',
    'content/Food4.png',
    'content/Sports1.png',
    'content/Sports2.png',
    'content/Sports3.png',
    'content/Sports4.png',
    'content/Sports5.png',
    'content/Sports6.png',
    'content/Barbecues1.png',
    'content/Barbecues2.png',
    'content/Barbecues3.png',
    'content/goofy1.png',
    'content/goofy2.png',
    'content/goofy3.png',
    'content/goofy4.png',
    'content/goofy5.png',
    'content/goofy6.png',
  ];

  // Preload images in batches to avoid overwhelming the browser
  function preloadImages() {
    const preloadBatchSize = 5;
    let currentBatch = 0;

    function loadBatch() {
      const start = currentBatch * preloadBatchSize;
      const end = Math.min(start + preloadBatchSize, imagesToPreload.length);

      for (let i = start; i < end; i++) {
        const img = new Image();
        img.src = imagesToPreload[i];
      }

      currentBatch++;

      // Load next batch after a short delay
      if (end < imagesToPreload.length) {
        setTimeout(loadBatch, 100);
      } else {
        console.log('All images preloaded');
      }
    }

    loadBatch();
  }

  // Start preloading when page is fully loaded
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', preloadImages);
  } else {
    preloadImages();
  }

  // Also preload on window load to catch any lazy images
  window.addEventListener('load', preloadImages);
})();

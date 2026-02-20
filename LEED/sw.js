// Service Worker for caching images and assets
const CACHE_NAME = 'bad-lab-v1';
const URLS_TO_CACHE = [
  '/',
  '/index.html',
  '/research.html',
  '/team.html',
  '/activities.html',
  '/publications.html',
  '/contact.html',
  '/styles.css',
  '/fetch-publications.js',
  // Team images
  '/content/team/salahuddin.jpg',
  '/content/team/KoushikDas.jpg',
  '/content/team/ChiragGarg.jpg',
  '/content/team/RamitDutta.jpg',
  '/content/team/RichardMcManus.png',
  '/content/team/KimPham.jpg',
  '/content/team/Wei-YangWeng.jpg',
  '/content/team/MyeongseopSong.jpg',
  '/content/team/JunghyeonHwang.jpg',
  // Research icons
  '/content/icon-materials.png',
  '/content/icon-simulation.png',
  '/content/icon-engineering.png',
  '/content/icon-engineering2.png',
  '/content/icon-computing.png',
  // Activity images
  '/content/Food1.png',
  '/content/Food2.png',
  '/content/Food3.png',
  '/content/Food4.png',
  '/content/Sports1.png',
  '/content/Sports2.png',
  '/content/Sports3.png',
  '/content/Sports4.png',
  '/content/Sports5.png',
  '/content/Sports6.png',
  '/content/Barbecues1.png',
  '/content/Barbecues2.png',
  '/content/goofy1.png',
  '/content/goofy2.png',
  '/content/goofy3.png',
  '/content/goofy4.png',
  '/content/goofy5.png',
  '/content/goofy6.png',
];

// Install event - cache all specified URLs
self.addEventListener('install', (event) => {
  event.waitUntil(
    caches.open(CACHE_NAME)
      .then((cache) => {
        console.log('Caching BAD Lab resources');
        return cache.addAll(URLS_TO_CACHE).catch((err) => {
          console.warn('Some resources could not be cached:', err);
          // Continue even if some resources fail to cache
          return Promise.resolve();
        });
      })
      .then(() => self.skipWaiting())
  );
});

// Activate event - clean up old caches
self.addEventListener('activate', (event) => {
  event.waitUntil(
    caches.keys().then((cacheNames) => {
      return Promise.all(
        cacheNames.map((cacheName) => {
          if (cacheName !== CACHE_NAME) {
            console.log('Deleting old cache:', cacheName);
            return caches.delete(cacheName);
          }
        })
      );
    }).then(() => self.clients.claim())
  );
});

// Fetch event - serve from cache, fallback to network
self.addEventListener('fetch', (event) => {
  const { request } = event;

  // Only handle GET requests
  if (request.method !== 'GET') {
    return;
  }

  // Skip non-GET requests and cross-origin requests
  if (request.mode === 'navigate') {
    event.respondWith(
      caches.match(request)
        .then((response) => response || fetch(request))
        .catch(() => new Response('Offline'))
    );
    return;
  }

  // For all other requests, try cache first, then network
  event.respondWith(
    caches.match(request)
      .then((response) => {
        if (response) {
          return response;
        }

        return fetch(request).then((response) => {
          // Don't cache non-successful responses
          if (!response || response.status !== 200 || response.type === 'error') {
            return response;
          }

          // Clone the response and cache it
          const responseToCache = response.clone();
          caches.open(CACHE_NAME)
            .then((cache) => {
              cache.put(request, responseToCache);
            });

          return response;
        });
      })
      .catch(() => {
        // Fallback for failed requests
        return caches.match(request);
      })
  );
});

// Fetch recent publications from Semantic Scholar with retry logic
async function fetchPublications() {
    const container = document.getElementById('recent-publications');
    if (!container) return;

    const maxRetries = 20;
    let attempt = 0;

    async function attemptFetch() {
        try {
            attempt++;
            // Search for Sayeef Salahuddin on Semantic Scholar
            const searchResponse = await fetch(
                'https://api.semanticscholar.org/graph/v1/author/search?query=Sayeef%20Salahuddin&limit=1'
            );
            const searchData = await searchResponse.json();

            if (!searchData.data || searchData.data.length === 0) {
                console.error('Author not found');
                if (attempt < maxRetries) {
                    console.log(`Retry attempt ${attempt}/${maxRetries} in 2 seconds...`);
                    await new Promise(resolve => setTimeout(resolve, 2000));
                    return attemptFetch();
                }
                return;
            }

            const authorId = searchData.data[0].authorId;

            // Fetch author's papers - request more papers to ensure we get recent ones
            const papersResponse = await fetch(
                `https://api.semanticscholar.org/graph/v1/author/${authorId}/papers?limit=200&fields=paperId,title,year,authors,venue,url,externalIds`
            );
            const papersData = await papersResponse.json();

            if (!papersData.data || papersData.data.length === 0) {
                console.error('No papers found');
                if (attempt < maxRetries) {
                    console.log(`Retry attempt ${attempt}/${maxRetries} in 2 seconds...`);
                    await new Promise(resolve => setTimeout(resolve, 2000));
                    return attemptFetch();
                }
                return;
            }

            // Filter for papers from 2015 onwards, sort by year (most recent first) and take top 5
            const publications = papersData.data
                .filter(paper => (paper.year || 0) >= 2015)
                .sort((a, b) => (b.year || 0) - (a.year || 0))
                .slice(0, 5);

            // If no publications found after filtering, retry
            if (publications.length === 0) {
                console.error('No publications found after filtering (attempt ' + attempt + '/' + maxRetries + ')');
                if (attempt < maxRetries) {
                    console.log(`Retry attempt ${attempt}/${maxRetries} in 2 seconds...`);
                    await new Promise(resolve => setTimeout(resolve, 2000));
                    return attemptFetch();
                }
            }

            // Clear existing placeholder content
            container.innerHTML = '';

            // Create publication items
            publications.forEach(pub => {
                const pubItem = document.createElement('div');
                pubItem.className = 'pub-item fade-in';

                const authors = pub.authors && pub.authors.length > 0
                    ? pub.authors.map(a => a.name).join(', ')
                    : 'Unknown authors';

                const venue = pub.venue || 'Publication';
                const year = pub.year || 'N/A';
                const googleScholarUrl = pub.externalIds?.CorpusId 
                    ? `https://scholar.google.com/scholar?q=${encodeURIComponent(pub.title)}`
                    : null;

                let html = `
                    <h3>${pub.title}</h3>
                    <p>${authors}</p>
                    <p class="pub-venue">${venue} (${year})</p>
                `;

                if (pub.url || googleScholarUrl) {
                    html += `<a href="${pub.url || googleScholarUrl}" target="_blank" class="pub-link">Read Paper →</a>`;
                }

                pubItem.innerHTML = html;
                container.appendChild(pubItem);
            });

            // Add "View All" link
            const viewAllItem = document.createElement('div');
            viewAllItem.className = 'pub-item view-all fade-in';
            viewAllItem.innerHTML = `<a href="https://scholar.google.com/citations?hl=en&user=Ay4FZFUAAAAJ&view_op=list_works&sortby=pubdate" class="pub-link" target="_blank">View complete publication list on Google Scholar →</a>`;
            container.appendChild(viewAllItem);

        } catch (error) {
            console.error(`Error fetching publications (attempt ${attempt}/${maxRetries}):`, error);
            if (attempt < maxRetries) {
                console.log(`Retrying in 2 seconds...`);
                await new Promise(resolve => setTimeout(resolve, 2000));
                return attemptFetch();
            } else {
                // Fallback to showing an error message after all retries failed
                container.innerHTML = '<p>Unable to load publications. Please try again later.</p>';
            }
        }
    }

    return attemptFetch();
}

// Load publications when DOM is ready
document.addEventListener('DOMContentLoaded', fetchPublications);

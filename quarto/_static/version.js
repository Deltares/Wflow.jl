function checkPathExists(url) {
    return new Promise((resolve, reject) => {
        var xhr = new XMLHttpRequest();
        xhr.open('HEAD', url, true);
        xhr.onreadystatechange = function() {
            if (xhr.readyState === 4) {
                if (xhr.status === 200) {
                    resolve(true);
                } else if (xhr.status === 404) {
                    resolve(false);
                } else {
                    reject(new Error(xhr.statusText));
                }
            }
        };
        xhr.onerror = function() {
            reject(new Error('Network Error'));
        };
        xhr.send();
    });
}

window.onload = function() {
    // Assuming you have a ul element in your HTML like this:
    // <ul id="myDropdown"></ul>

    // Fetch the JSON data
    fetch("https://raw.githubusercontent.com/Deltares/Delft-FIAT/gh-pages/switcher.json")
      .then(response => response.json())
      .then(data => {
        console.log('Data loaded:', data); // Log the loaded data

        const dropdown = document.querySelector('#nav-menu-version').nextElementSibling;
        console.log('Dropdown element:', dropdown); // Log the dropdown element

        // Clear all existing dropdown items
        dropdown.innerHTML = '';

        data.forEach(item => {
            console.log('Adding item:', item); // Log the item being added

            // Create a new li element
            const li = document.createElement('li');

            // Create a new a element
            const a = document.createElement('a');
            a.className = 'dropdown-item';
            a.href = item.url; // Use the 'url' property as the href
            a.textContent = item.name; // Use the 'name' property as the text

            // Add the a element to the li
            li.appendChild(a);

            // Add the li to the dropdown
            dropdown.appendChild(li);
        });

        console.log('Dropdown after adding items:', dropdown); // Log the dropdown after adding items

        // Get all dropdown items within the specific dropdown menu
        var dropdownMenu = document.querySelector('#nav-menu-version').nextElementSibling;

        var dropdownItems = dropdownMenu.querySelectorAll('.dropdown-item');

        // Get the current page in chunks
        var currentPagePath = window.location.pathname.split('/');

        for (var i = 0; i < dropdownItems.length; i++) {
            // Get textcontent
            var textContent = dropdownItems[i].textContent;

            // Get the index of the current version
            var index = currentPagePath.indexOf(textContent);

            if (index !== -1) {
                // Remove the active-item class from all items
                for (var j = 0; j < dropdownItems.length; j++) {
                    dropdownItems[j].classList.remove('active-item');
                }

                dropdownItems[i].classList.add('active-item');
                break
            }
        }

        console.log('current page path', currentPagePath);

        // Loop through each dropdown item
        for (var i = 0; i < dropdownItems.length; i++) {
            // Add click event listener to each item
            dropdownItems[i].addEventListener('click', function(event) {
                // Prevent default action
                event.preventDefault();

                // Get the clicked item's text
                var itemText = this.textContent;
                // var itemHref = this.getAttribute('href')

                // Loop through each dropdown item again to find a match in the current page's path
                for (var j = 0; j < dropdownItems.length; j++) {
                    // Get the dropdown item's text
                    var dropdownText = dropdownItems[j].textContent;
                    console.log('Dropdown item:', dropdownText);

                    // Find the index of the dropdownText in the current page's path
                    var index = currentPagePath.indexOf(dropdownText);

                    // If the dropdownText is found in the current page's path
                    if (index !== -1) {
                        // Construct the new URL relative to the dropdownText and append the itemText
                        addElements = currentPagePath.slice(index + 1, )
                        relativePath = '../'.repeat(addElements.length)
                        var newUrl = relativePath + itemText + '/' + addElements.join('/')
                        console.log('Clicked item:', newUrl);

                        // Redirect to the new URL
                        checkPathExists(newUrl)
                            .then(exists => {
                                if (exists) {
                                    window.location.href = newUrl;
                                } else {
                                    console.log('Path does not exist, referring to home page');
                                    window.location.href = relativePath + itemText + '/';
                                }
                            })

                        // Exit the loop
                        break;
                    }
                }
            });
        }

        })
        .catch(error => console.error('Error:', error)); // Log any errors
}

import ApiError from '../shared/xhr/ApiError';
import fetch from '../shared/xhr/fetch';

async function checkValidity(smiles: string) {
  const response = await fetch('api/valid/', {
    method: 'POST',
    body: JSON.stringify({ smiles }),
  });

  if (!response.ok) {
    throw new ApiError(response.clone());
  }

  return response.json();
}

async function checkProps(smiles: string) {
  const response = await fetch('api/number/', {
    method: 'POST',
    body: JSON.stringify({ smiles }),
  });

  if (!response.ok) {
    throw new ApiError(response.clone());
  }

  return response.json();
}

async function getImage(smiles: string) {
  const response = await fetch('api/image/', {
    method: 'GET',
    body: JSON.stringify({ smiles }),
  });

  if (!response.ok) {
    throw new ApiError(response.clone());
  }

  return response.json();
}

declare global {
  interface Window {
    checkValidity: any;
    getImage: any;
    submit: any;
  }
}
window.checkValidity = checkValidity;
export { checkValidity };

// async function registerNewRGroup(smiles: string) {
//   const response = await fetch(rgroupApiUrls.create, {
//     method: 'POST',
//     body: JSON.stringify({ smiles }),
//   });

//   if (!response.ok) {
//     throw new ApiError(response.clone());
//   }

//   return response.json() as Promise<RGroup>;
// }

// async function enumerateProperties(enumerationOpts: {
//   core_smiles: string;
//   rgroup_smiles: Record<string, string[]>;
// }) {
//   const response = await fetch(enumerationApiUrls.create, {
//     method: 'POST',
//     body: JSON.stringify(enumerationOpts),
//   });

//   if (!response.ok) {
//     throw new ApiError(response.clone());
//   }

//   return response.json() as Promise<EnumerateResponse>;
// }

// export type { RGroup, EnumerateResponse };
// export { fetchAllRGroups, registerNewRGroup, enumerateProperties };

import {HttpClient} from '@angular/common/http';
import {Component, ViewChild, AfterViewInit, inject} from '@angular/core';
import {MatPaginator, MatPaginatorModule} from '@angular/material/paginator';
import {MatSort, MatSortModule, SortDirection} from '@angular/material/sort';
import {merge, Observable, of as observableOf} from 'rxjs';
import {catchError, map, startWith, switchMap} from 'rxjs/operators';
import {MatTableModule} from '@angular/material/table';
import {MatProgressSpinnerModule} from '@angular/material/progress-spinner';
import {MaveDBData, MaveDBFilters, MaveDBResponse} from "../model/mavedb";
import {MatFormFieldModule} from "@angular/material/form-field";
import {MatInputModule} from "@angular/material/input";
import {MatCardModule} from "@angular/material/card";
import {MatListModule} from "@angular/material/list";
import {MatDividerModule} from "@angular/material/divider";

@Component({
  selector: 'app-data-portal',
  imports: [
    MatProgressSpinnerModule,
    MatTableModule,
    MatSortModule,
    MatPaginatorModule,
    MatFormFieldModule,
    MatInputModule,
    MatCardModule,
    MatListModule,
    MatDividerModule
  ],
  templateUrl: './data-portal.component.html',
  styleUrl: './data-portal.component.scss'
})
export class DataPortalComponent implements AfterViewInit {

  private _httpClient = inject(HttpClient);

  displayedColumns: string[] = [
    'urn',
    'sequenceType',
    'geneName',
    'geneCategory',
    'publicationYear',
    'numVariants'
  ];
  exampleDatabase: MaveDBDataService | null | undefined;
  data: MaveDBData[] = [];
  aggregations: any;
  searchValue!: string;
  filters: MaveDBFilters = {
    sequenceType: [],
    geneCategory: [],
    publicationYear: []
  };

  resultsLength = 0;
  isLoadingResults = true;

  @ViewChild(MatPaginator) paginator!: MatPaginator;
  @ViewChild(MatSort) sort!: MatSort;

  ngAfterViewInit() {
    this.exampleDatabase = new MaveDBDataService(this._httpClient);

    // If the user changes the sort order, reset back to the first page.
    this.sort.sortChange.subscribe(() => (this.paginator.pageIndex = 0));

    merge(this.sort.sortChange, this.paginator.page)
      .pipe(
        startWith({}),
        switchMap(() => {
          this.isLoadingResults = true;
          return this.exampleDatabase!.getAllMaveDBData(
            this.paginator.pageIndex,
            this.paginator.pageSize,
            this.sort.active,
            this.sort.direction,
            this.searchValue,
            this.filters,
          ).pipe(catchError(() => observableOf(null)));
        }),
        map(data => {
          // Flip flag to show that loading has finished.
          this.isLoadingResults = false;

          if (data === null) {
            return [];
          }

          // Only refresh the result length if there is new data. In case of rate
          // limit errors, we do not want to reset the paginator to zero, as that
          // would prevent users from re-triggering requests.
          this.resultsLength = data.total;
          this.aggregations = data.aggregations;
          return data.results;
        }),
      )
      .subscribe(data => (this.data = data));
  }

  search(event: Event) {
    this.searchValue = (event.target as HTMLInputElement).value.trim().toLowerCase();
    this.paginator.page.emit();
  }

  applyFilter(filterKey: string, filterValue: string) {
    // @ts-ignore
    const index = this.filters[filterKey].indexOf(filterValue);
    if (index > -1) {
      // @ts-ignore
      this.filters[filterKey].splice(index, 1)
    } else {
      // @ts-ignore
      this.filters[filterKey].push(filterValue);
    }
    this.paginator.page.emit();
  }

}

/** An example database that the data source uses to retrieve data for the table. */
export class MaveDBDataService {
  constructor(private _httpClient: HttpClient) {}

  getAllMaveDBData(
    start: number,
    size: number,
    sortField: string = 'publicationYear',
    sortOrder: string = 'desc',
    query: any = null,
    filters: MaveDBFilters): Observable<MaveDBResponse> {
    const requestUrl = 'https://perturbation-catalogue-be-959149465821.europe-west2.run.app/search';
    let params = {
      start: start*size,
      size: size,
      sort_field: sortField,
      sort_order: sortOrder};

    if (query) {
      // @ts-ignore
      params['q'] = query;
    }

    let include_filters = false;
    let filter = '';
    if (filters.publicationYear.length > 0) {
      include_filters = true;
      filter = `publicationYear:${filters.publicationYear}`;
    }
    if (filters.geneCategory.length > 0) {
      if (include_filters) {
        filter = `${filter}+geneCategory:${filters.geneCategory}`;
      } else {
        include_filters = true;
        filter = `geneCategory:${filters.geneCategory}`;
      }
    }
    if (filters.sequenceType.length > 0) {
      if (include_filters) {
        filter = `${filter}+sequenceType:${filters.sequenceType}`;
      } else {
        include_filters = true;
        filter = `sequenceType:${filters.sequenceType}`;
      }
    }
    if (include_filters) {
      // @ts-ignore
      params['data_filter'] = filter;
    }
    return this._httpClient.get<MaveDBResponse>(requestUrl, {params: params});
  }
}
